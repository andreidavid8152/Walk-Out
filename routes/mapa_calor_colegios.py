# =========================================================
# 1. IMPORTACIONES Y CONFIGURACIÓN BÁSICA
# =========================================================
from flask import Blueprint, render_template, request
import geopandas as gpd
import pandas as pd
import folium
from folium.plugins.treelayercontrol import TreeLayerControl
from shapely.geometry import box
import os
from shapely.geometry import Point
from branca.colormap import linear

mapa_colegios_bp = Blueprint("mapa_calor_colegios", __name__)

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, "..", "data")

# Rutas de archivos
EXCEL_PATH = os.path.join(DATA_DIR, "universidades_colegios.xlsx")
SHEET_COL = "Colegios"
CSV_EST = os.path.join(DATA_DIR, "ubicacionEstudiantesPeriodo.csv")
GJSON_RURAL = os.path.join(DATA_DIR, "parroquiasRurales.geojson")
GJSON_URB = os.path.join(DATA_DIR, "parroquiasUrbanas.geojson")
GJSON_BUSES = os.path.join(DATA_DIR, "estacionesBuses.geojson")
GJSON_METRO = os.path.join(DATA_DIR, "estacionesMetro.geojson")
GJSON_PARADAS = os.path.join(DATA_DIR, "paradasBuses.geojson")


# =========================================================
# 2. RUTA PRINCIPAL DEL MAPA
# =========================================================
@mapa_colegios_bp.route("/mapacalor/colegios")
def mapa():
    # -----------------------------------------------------------------
    # 2-A. Parámetro de período (solo para mantener tu selector)
    # -----------------------------------------------------------------
    selected_periodo = request.args.get("periodo")
    df_all = pd.read_csv(CSV_EST, sep=";").rename(columns={"Semestre": "periodo"})
    df_all["periodo"] = df_all["periodo"].astype(str)
    periodos = sorted(df_all["periodo"].unique())
    if selected_periodo not in periodos:
        selected_periodo = periodos[0]

    # -----------------------------------------------------------------
    # 2-B. Carga de datos geoespaciales
    # -----------------------------------------------------------------
    # Parroquias
    gdf_rurales = gpd.read_file(GJSON_RURAL).rename(columns={"DPA_DESPAR": "nombre"})
    gdf_urbanas = gpd.read_file(GJSON_URB).rename(columns={"dpa_despar": "nombre"})
    gdf_rurales["tipo"] = "rural"
    gdf_urbanas["tipo"] = "urbana"
    gdf_parroquias = pd.concat(
        [
            gdf_rurales[["nombre", "geometry", "tipo"]],
            gdf_urbanas[["nombre", "geometry", "tipo"]],
        ],
        ignore_index=True,
    ).set_crs("EPSG:4326")

    # Transporte
    gdf_buses = gpd.read_file(GJSON_BUSES).to_crs("EPSG:4326")
    gdf_metro = gpd.read_file(GJSON_METRO).to_crs("EPSG:4326")
    gdf_paradas = gpd.read_file(GJSON_PARADAS).to_crs("EPSG:4326")

    # Colegios AAA
    df_col = pd.read_excel(EXCEL_PATH, sheet_name=SHEET_COL).rename(columns=lambda c: c.strip())
    df_col["LATITUD"] = pd.to_numeric(df_col["LATITUD"], errors="coerce")
    df_col["LONGITUD"] = pd.to_numeric(df_col["LONGITUD"], errors="coerce")
    df_col.dropna(subset=["LATITUD", "LONGITUD"], inplace=True)

    df_aaa = df_col[df_col["TIPO"].str.upper() == "AAA"]

    # Crear GeoDataFrame de colegios AAA
    gdf_colegios = gpd.GeoDataFrame(
        df_aaa,
        geometry=gpd.points_from_xy(df_aaa["LONGITUD"], df_aaa["LATITUD"]),
        crs="EPSG:4326"
    )

    # -----------------------------------------------------------------
    # 2-C. Preparación de la grilla basada en la mediana del área de parroquias
    # -----------------------------------------------------------------

    # Reproyectar parroquias a metros para calcular áreas
    gdf_parroquias_m = gdf_parroquias.to_crs(epsg=32717)
    gdf_parroquias_m["area_m2"] = gdf_parroquias_m.geometry.area

    # Calcular mediana de áreas
    mediana_area_m2 = gdf_parroquias_m["area_m2"].median()
    lado_celda_m = mediana_area_m2**0.5  # Lado de celda cuadrada en metros

    # Crear grilla en CRS métrico
    minx, miny, maxx, maxy = gdf_parroquias_m.total_bounds

    grid_cells = []
    x = minx
    while x < maxx:
        y = miny
        while y < maxy:
            cell = box(x, y, x + lado_celda_m, y + lado_celda_m)
            grid_cells.append(cell)
            y += lado_celda_m
        x += lado_celda_m

    gdf_grilla = gpd.GeoDataFrame(geometry=grid_cells, crs="EPSG:32717")

    # Volver a EPSG:4326 para el mapa y futuros joins
    gdf_grilla = gdf_grilla.to_crs("EPSG:4326")

    # ─── 2-D. CÁLCULO DE DENSIDAD EN LA GRILLA ───────────────────────────
    # 1) Crear GeoDataFrames de puntos

    gdf_buses_points = gdf_buses.copy()
    gdf_buses_points = gdf_buses_points.to_crs("EPSG:32717")
    gdf_buses_points.geometry = gdf_buses_points.geometry.centroid
    gdf_buses_points = gdf_buses_points.set_geometry("geometry").to_crs("EPSG:4326")

    gdf_metro_points = gdf_metro.copy()
    gdf_metro_points = gdf_metro_points.to_crs("EPSG:32717")
    gdf_metro_points.geometry = gdf_metro_points.geometry.centroid
    gdf_metro_points = gdf_metro_points.set_geometry("geometry").to_crs("EPSG:4326")

    # las paradas ya son puntos
    gdf_paradas_points = gdf_paradas.copy()

    # 2) Unir todos los puntos en un solo GeoDataFrame
    gdf_puntos = pd.concat(
        [
            gdf_buses_points[["geometry"]],
            gdf_metro_points[["geometry"]],
            gdf_paradas_points[["geometry"]],
            gdf_colegios[["geometry"]],
        ],
        ignore_index=True,
    ).set_crs("EPSG:4326")

    # 3) Spatial join para asignar cada punto a su celda
    join = gpd.sjoin(gdf_grilla, gdf_puntos, how="left", predicate="contains")

    join_valid = join[join["index_right"].notnull()]

    # 4) Contar puntos por celda y guardar en una nueva columna
    counts = join_valid.groupby(join_valid.index).size()
    gdf_grilla["count"] = counts.reindex(gdf_grilla.index).fillna(0).astype(int)

    # 5) Crear un colormap de YlOrRd según el rango de 'count'
    colormap = linear.YlOrRd_09.scale(0, gdf_grilla["count"].max())
    colormap.caption = "Densidad de puntos de interés"
    # ────────────────────────────────────────────────────────────────

    # ================================================================
    # 3. CONSTRUCCIÓN DEL MAPA Y CAPAS
    # ================================================================
    m = folium.Map(location=[-0.20, -78.50], zoom_start=11, tiles="cartodbpositron")

    # ─── 3-0. Capa de Mapa de Calor (cloropleth sobre la grilla) ─────────────────
    fg_heat = folium.FeatureGroup(name="Mapa de Calor", show=True).add_to(m)
    folium.GeoJson(
        gdf_grilla,
        style_function=lambda feature: {
            "fillColor": (
                colormap(feature["properties"]["count"])
                if feature["properties"]["count"] > 0
                else "white"
            ),
            "color": "grey",
            "weight": 0.6,
            "fillOpacity": 0.7,
        },
        tooltip=folium.GeoJsonTooltip(
            fields=["count"], aliases=["Total puntos:"], localize=True
        ),
    ).add_to(fg_heat)

    # Añadir la leyenda (colormap) al mapa
    colormap.add_to(m)
    # ─────────────────────────────────────────────────────────────────────────

    # --- 3-A. Parroquias -------------------------------------------
    fg_parroquias = folium.FeatureGroup(name="Parroquias").add_to(m)
    for _, row in gdf_parroquias.iterrows():
        folium.GeoJson(
            {
                "type": "Feature",
                "geometry": row.geometry.__geo_interface__,
                "properties": {"nombre": row["nombre"]},
            },
            style_function=lambda _: {
                "fillColor": "white",
                "color": "black",
                "weight": 1,
                "fillOpacity": 0.01,
            },
            tooltip=folium.GeoJsonTooltip(fields=["nombre"], aliases=["Parroquia:"]),
        ).add_to(fg_parroquias)

    # --- 3-B. Colegios AAA -------------------------------------------
    fg_colegios_aaa = folium.FeatureGroup(name="Colegios AAA").add_to(m)

    for _, row in df_aaa.iterrows():
        folium.Marker(
            location=[row["LATITUD"], row["LONGITUD"]],
            tooltip=row["COLEGIO"].strip(),
            icon=folium.Icon(color="blue", icon="graduation-cap", prefix="fa"),
        ).add_to(fg_colegios_aaa)

    # --- 3-C. Transporte público ------------------------------------
    ## Estaciones de buses
    fg_buses = folium.FeatureGroup(name="Estaciones de Buses", show=False).add_to(m)
    for _, row in gdf_buses.iterrows():
        geom = row.geometry
        folium.GeoJson(
            geom.__geo_interface__,
            style_function=lambda _: {
                "fillColor": "red",
                "color": "red",
                "weight": 1.5,
                "fillOpacity": 0.4,
            },
            tooltip="Estación de Bus",
        ).add_to(fg_buses)
        folium.Marker(
            location=[geom.centroid.y, geom.centroid.x],
            icon=folium.Icon(icon="bus", prefix="fa", color="red"),
            tooltip="Estación de Bus",
        ).add_to(fg_buses)

    ## Estaciones de metro
    fg_metro = folium.FeatureGroup(name="Estaciones de Metro", show=False).add_to(m)
    for _, row in gdf_metro.iterrows():
        geom = row.geometry
        nombre_estacion = f"Estación de metro: {row.get('nam', 'Desconocida')}"
        folium.GeoJson(
            geom.__geo_interface__,
            style_function=lambda _: {
                "fillColor": "purple",
                "color": "purple",
                "weight": 1.5,
                "fillOpacity": 0.4,
            },
            tooltip=nombre_estacion,
        ).add_to(fg_metro)
        folium.Marker(
            location=[geom.centroid.y, geom.centroid.x],
            icon=folium.Icon(icon="subway", prefix="fa", color="purple"),
            tooltip=nombre_estacion,
        ).add_to(fg_metro)

    ## Paradas de buses (puntos)
    fg_paradas = folium.FeatureGroup(name="Paradas de Buses", show=False).add_to(m)
    for _, row in gdf_paradas.iterrows():
        punto = row.geometry
        folium.CircleMarker(
            location=[punto.y, punto.x],
            radius=4,
            color="darkgreen",
            fill=True,
            fill_color="limegreen",
            fill_opacity=0.8,
            tooltip="Parada de Bus",
        ).add_to(fg_paradas)

    # ================================================================
    # 4. CONTROL DE CAPAS (TreeLayerControl)
    # ================================================================
    TreeLayerControl(
        overlay_tree=[
            {"label": "Parroquias", "layer": fg_parroquias},
            {"label": "Colegios AAA", "layer": fg_colegios_aaa},
            {
                "label": "Transporte Público",
                "select_all_checkbox": True,
                "children": [
                    {"label": "Estaciones de Buses", "layer": fg_buses},
                    {"label": "Estaciones de Metro", "layer": fg_metro},
                    {"label": "Paradas de Buses", "layer": fg_paradas},
                ],
            },
            {"label": "Mapa de Calor", "layer": fg_heat},
        ]
    ).add_to(m)

    # ================================================================
    # 5. RENDERIZACIÓN DE LA PLANTILLA
    # ================================================================
    return render_template(
        "mapa_calor_colegios.html",
        mapa=m.get_root().render(),
        map_name=m.get_name(),
        periodos=periodos,
        selected_periodo=selected_periodo,
        ruta_activa="colegios"  
    )
