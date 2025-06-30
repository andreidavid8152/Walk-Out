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

mapa_hospitales_bp = Blueprint("mapa_calor_hospitales", __name__)

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, "..", "data")

# Rutas de archivos - Solo los necesarios
GJSON_RURAL       = os.path.join(DATA_DIR, "parroquiasRurales.geojson")
GJSON_URB         = os.path.join(DATA_DIR, "parroquiasUrbanas.geojson")
GJSON_PROVINCIAS = os.path.join(DATA_DIR, "ec.json")

# Hospitales desde Excel
HOSPITALES_XLSX = os.path.join(DATA_DIR, "Hospitales.xlsx")

# =========================================================
# 2. RUTA PRINCIPAL DEL MAPA
# =========================================================
@mapa_hospitales_bp.route("/mapacalor/hospitales")
def mapa():
    # -----------------------------------------------------------------
    # 2-A. Parámetro de período (solo para mantener tu selector)
    # -----------------------------------------------------------------
    selected_periodo = request.args.get("periodo")
    # Como no tenemos CSV de estudiantes, usamos un período por defecto
    periodos = ["202410", "202420"]
    if selected_periodo not in periodos:
        selected_periodo = periodos[0]

    # -----------------------------------------------------------------
    # 2-B. Carga de datos geoespaciales
    # -----------------------------------------------------------------
    # Provincias
    gdf_provincias = gpd.read_file(GJSON_PROVINCIAS)
    gdf_provincias = gdf_provincias[gdf_provincias["name"].str.lower() == "pichincha"]
    print(gdf_provincias.columns)
    print(gdf_provincias.head())

    # Hospitales
    df_hosp = pd.read_excel(HOSPITALES_XLSX)
    df_hosp.columns = df_hosp.columns.str.strip().str.upper()
    gdf_hospitales = gpd.GeoDataFrame(
        df_hosp,
        geometry=gpd.points_from_xy(df_hosp["LONGITUD"], df_hosp["LATITUD"]),
        crs="EPSG:4326"
    )

    # ================================================================
    # 3. CONSTRUCCIÓN DEL MAPA Y CAPAS
    # ================================================================
    m = folium.Map(location=[-0.1807, -78.4678], zoom_start=10, tiles="cartodbpositron")

    # --- 3-A. Provincias -------------------------------------------
    fg_provincias = folium.FeatureGroup(name="Provincias").add_to(m)
    folium.GeoJson(
        gdf_provincias,
        style_function=lambda _: {"fillColor": "white", "color": "black", "weight": 1, "fillOpacity": 0.01},
        tooltip=folium.GeoJsonTooltip(fields=["name"], aliases=["Provincia:"]),
    ).add_to(fg_provincias)

    # --- 3-B. Hospitales y Centros de Salud con filtros avanzados ---
    # Definir colores
    COLOR_CONVENIO = "green" 
    COLOR_SIN_CONVENIO = "red"

    # Filtrar hospitales y centros de salud
    tipos = gdf_hospitales["TIPO"].str.upper().unique()
    institucion_col = "INSTITUCION"

    # Crear estructura de capas
    overlay_tree = [
        {"label": "Pichincha", "layer": fg_provincias},
    ]

    for tipo in ["HOSPITAL", "CENTRO DE SALUD"]:
        fg_tipo = folium.FeatureGroup(name=tipo.title(), show=(tipo=="HOSPITAL")).add_to(m)
        hijos = []
        # Públicos
        fg_publico = folium.FeatureGroup(name=f"{tipo.title()}s Públicos").add_to(m)
        # --- Con Convenio Públicos ---
        fg_publico_convenio = folium.FeatureGroup(name=f"{tipo.title()}s Públicos Con Convenio").add_to(m)
        filtro_convenio_publico = (
            (gdf_hospitales["TIPO"].str.upper() == tipo)
            & (gdf_hospitales["INSTITUCION"].str.upper() == "PUBLICO")
            & (gdf_hospitales["CONVENIO"].str.upper() == "SI")
        )
        for _, row in gdf_hospitales[filtro_convenio_publico].iterrows():
            lat, lon = row.geometry.y, row.geometry.x
            nombre = row.get("NOMBRE", tipo.title())
            folium.Marker(
                location=[lat, lon],
                icon=folium.Icon(icon="hospital", prefix="fa", color=COLOR_CONVENIO),
                tooltip=nombre,
            ).add_to(fg_publico_convenio)
        fg_publico.add_child(fg_publico_convenio)
        # --- Sin Convenio Públicos ---
        fg_publico_sinconvenio = folium.FeatureGroup(name=f"{tipo.title()}s Públicos Sin Convenio").add_to(m)
        filtro_sinconvenio_publico = (
            (gdf_hospitales["TIPO"].str.upper() == tipo)
            & (gdf_hospitales["INSTITUCION"].str.upper() == "PUBLICO")
            & (gdf_hospitales["CONVENIO"].str.upper() == "NO")
        )
        for _, row in gdf_hospitales[filtro_sinconvenio_publico].iterrows():
            lat, lon = row.geometry.y, row.geometry.x
            nombre = row.get("NOMBRE", tipo.title())
            folium.Marker(
                location=[lat, lon],
                icon=folium.Icon(icon="hospital", prefix="fa", color=COLOR_SIN_CONVENIO),
                tooltip=nombre,
            ).add_to(fg_publico_sinconvenio)
        fg_publico.add_child(fg_publico_sinconvenio)
        hijos.append({
            "label": "Públicos",
            "layer": fg_publico,
            "children": [
                {"label": "Con Convenio", "layer": fg_publico_convenio},
                {"label": "Sin Convenio", "layer": fg_publico_sinconvenio},
            ]
        })
        # Privado
        fg_privado = folium.FeatureGroup(name=f"{tipo.title()}s Privados").add_to(m)
        # --- Con Convenio Privados ---
        fg_privado_convenio = folium.FeatureGroup(name=f"{tipo.title()}s Privados Con Convenio").add_to(m)
        filtro_convenio_privado = (
            (gdf_hospitales["TIPO"].str.upper() == tipo)
            & (gdf_hospitales["INSTITUCION"].str.upper() == "PRIVADA")
            & (gdf_hospitales["CONVENIO"].str.upper() == "SI")
        )
        for _, row in gdf_hospitales[filtro_convenio_privado].iterrows():
            lat, lon = row.geometry.y, row.geometry.x
            nombre = row.get("NOMBRE", tipo.title())
            folium.Marker(
                location=[lat, lon],
                icon=folium.Icon(icon="hospital", prefix="fa", color=COLOR_CONVENIO),
                tooltip=nombre,
            ).add_to(fg_privado_convenio)
        fg_privado.add_child(fg_privado_convenio)
        # --- Sin Convenio Privados ---
        fg_privado_sinconvenio = folium.FeatureGroup(name=f"{tipo.title()}s Privados Sin Convenio").add_to(m)
        filtro_sinconvenio_privado = (
            (gdf_hospitales["TIPO"].str.upper() == tipo)
            & (gdf_hospitales["INSTITUCION"].str.upper() == "PRIVADA")
            & (gdf_hospitales["CONVENIO"].str.upper() == "NO")
        )
        for _, row in gdf_hospitales[filtro_sinconvenio_privado].iterrows():
            lat, lon = row.geometry.y, row.geometry.x
            nombre = row.get("NOMBRE", tipo.title())
            folium.Marker(
                location=[lat, lon],
                icon=folium.Icon(icon="hospital", prefix="fa", color=COLOR_SIN_CONVENIO),
                tooltip=nombre,
            ).add_to(fg_privado_sinconvenio)
        fg_privado.add_child(fg_privado_sinconvenio)
        hijos.append({
            "label": "Privados",
            "layer": fg_privado,
            "children": [
                {"label": "Con Convenio", "layer": fg_privado_convenio},
                {"label": "Sin Convenio", "layer": fg_privado_sinconvenio},
            ]
        })
        overlay_tree.append({
            "label": tipo.title(),
            "select_all_checkbox": True,
            "children": hijos,
        })

    # ================================================================
    # 4. CONTROL DE CAPAS (TreeLayerControl)
    # ================================================================
    TreeLayerControl(
        overlay_tree=overlay_tree
    ).add_to(m)

    # ================================================================
    # 5. RENDERIZACIÓN DE LA PLANTILLA
    # ================================================================
    return render_template(
        "mapa_calor_hospitales.html",
        mapa=m.get_root().render(),
        map_name=m.get_name(),
        periodos=periodos,
        selected_periodo=selected_periodo,
        ruta_activa="hospitales",
    )
