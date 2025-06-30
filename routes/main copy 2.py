from flask import Blueprint, render_template, request
import geopandas as gpd
import pandas as pd
import folium
from folium.plugins.treelayercontrol import TreeLayerControl
from shapely.geometry import box, Point
from shapely.strtree import STRtree
import random
import os
import json
import collections
from math import floor
from scipy.stats import chi2
import numpy as np
from shapely.affinity import translate
from shapely.ops import transform
import pyproj
from copy import deepcopy

main_bp = Blueprint("main", __name__)

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, "..", "data")

EXCEL_PATH = os.path.join(DATA_DIR, "universidades_colegios.xlsx")
SHEET_UNI = "Universidades"
CSV_EST = os.path.join(DATA_DIR, "ubicacionEstudiantesPeriodo.csv")
GJSON_RURAL = os.path.join(DATA_DIR, "parroquiasRurales.geojson")
GJSON_URB = os.path.join(DATA_DIR, "parroquiasUrbanas.geojson")
GJSON_BUSES = os.path.join(DATA_DIR, "estacionesBuses.geojson")
GJSON_METRO = os.path.join(DATA_DIR, "estacionesMetro.geojson")
GJSON_PARADAS_BUSES = os.path.join(DATA_DIR, "paradasBuses.geojson")
GJSON_ALIMENTADORES = os.path.join(DATA_DIR, "alimentadores.geojson")

# Cargar nombres de alimentadores desde JSON
with open(os.path.join(DATA_DIR, "idAlimentadores.json"), encoding="utf-8") as f:
    alimentador_nombre_map = {
        item["code"]: item["name"] for item in json.load(f)["codedValues"]
    }

@main_bp.route("/")
def mapa():
    # Fijar la semilla para reproducibilidad
    random.seed(42)
    np.random.seed(42)

    # 🎯 Incluye solo lógica de parroquias + universidades + filtros

    # Periodo
    selected_periodo = request.args.get("periodo")
    df_all = pd.read_csv(CSV_EST, sep=";").rename(columns={"Semestre": "periodo"})
    df_all["periodo"] = df_all["periodo"].astype(str)
    periodos = sorted(df_all["periodo"].unique())
    if selected_periodo not in periodos:
        selected_periodo = periodos[0]

    # Alimentadores
    gdf_alimentadores = gpd.read_file(GJSON_ALIMENTADORES).to_crs("EPSG:4326")

    # ============================================================
    # 🔹 1. CÁLCULO DE GRILLA PARA ALIMENTADORES
    # ============================================================

    gdf_alimentadores_m = gdf_alimentadores.to_crs(epsg=32717)
    gdf_alimentadores_m["area_m2"] = gdf_alimentadores_m.geometry.area


    # Dimensiones del área
    minx_a, miny_a, maxx_a, maxy_a = gdf_alimentadores_m.total_bounds
    ancho_total = maxx_a - minx_a
    alto_total = maxy_a - miny_a

    # 🔧 Número deseado de celdas (puedes cambiarlo a 2, 4, 16, etc.)
    NUM_CELDAS_DESEADAS = 7

    # 🔢 Calcula el lado para dividir el área total en ese número de celdas cuadradas
    area_total = ancho_total * alto_total
    lado_celda_alim = (area_total / NUM_CELDAS_DESEADAS) ** 0.5
    print(f"[🔍] Lado ajustado para dividir en {NUM_CELDAS_DESEADAS} celdas: {lado_celda_alim:.2f} m")

    filas = floor((maxy_a - miny_a) / lado_celda_alim)
    columnas = floor((maxx_a - minx_a) / lado_celda_alim)

    grid_cells_a = []
    for i in range(columnas):
        for j in range(filas):
            x0 = minx_a + i * lado_celda_alim
            y0 = miny_a + j * lado_celda_alim
            grid_cells_a.append(box(x0, y0, x0 + lado_celda_alim, y0 + lado_celda_alim))

    print(f"✅ Grilla generada: {filas} filas × {columnas} columnas = {len(grid_cells_a)} celdas")

    gdf_grilla_alimentadores = gpd.GeoDataFrame(
        geometry=grid_cells_a, crs="EPSG:32717"
    ).to_crs("EPSG:4326")

    # Intersección celda × alimentadores
    gdf_inter = gpd.overlay(
        gdf_grilla_alimentadores[["geometry"]],
        gdf_alimentadores[["geometry", "alimentadorid", "zonainfluencia"]],
        how="intersection",
    ).explode(ignore_index=True)

    gdf_inter = gdf_inter[gdf_inter.geometry.geom_type == "Polygon"]
    gdf_inter["area_m2"] = gdf_inter.to_crs("EPSG:32717").area
    gdf_inter = gdf_inter[gdf_inter["area_m2"] > 25]

    # Filtramos solo los que tienen un bloque válido
    gdf_inter = gdf_inter[gdf_inter["zonainfluencia"].notna()]

    # Reproyectamos para trabajar en metros
    gdf_inter_proj = gdf_inter.to_crs(epsg=32717)
    gdf_inter_proj["area_m2"] = gdf_inter_proj.geometry.area

    # Proyecta la grilla de alimentadores a metros y crea el STRtree
    gdf_grilla_alimentadores_proj = gdf_grilla_alimentadores.to_crs("EPSG:32717")
    tree = STRtree(gdf_grilla_alimentadores_proj.geometry)

    num_simulaciones = 1000
    best_chi2 = float('inf')
    best_centroides = None

    zonas = list(gdf_inter_proj.groupby("zonainfluencia"))

    def random_point_in_polygon(polygon):
        minx, miny, maxx, maxy = polygon.bounds
        for _ in range(30):  # Intenta hasta 30 veces
            p = Point(random.uniform(minx, maxx), random.uniform(miny, maxy))
            if polygon.contains(p):
                return p
        return polygon.representative_point()  # Fallback

    for _ in range(num_simulaciones):
        np.random.shuffle(zonas)  # Orden aleatorio de alimentadores
        conteo_celdas = {i: 0 for i in range(len(gdf_grilla_alimentadores_proj))}
        num_celdas = len(gdf_grilla_alimentadores_proj)
        num_centroides = len(zonas) * 2  # 2 centroides por alimentador
        esperada = num_centroides / num_celdas
        bloques_centroides = []
        for bloque, grupo in zonas:
            grupo_sorted = grupo.sort_values("area_m2", ascending=False)
            posibles = []
            for _, row in grupo_sorted.iterrows():
                centroide = row.geometry.representative_point()
                res = tree.query(centroide)
                if len(res) > 0:
                    celda_idx = res[0]
                    posibles.append((row, centroide, celda_idx))
            if len(posibles) < 2:
                posibles = posibles * 2
            # Ordena por celdas menos saturadas
            posibles = sorted(posibles, key=lambda x: conteo_celdas[x[2]])
            usados = set()
            centroides = []
            for row, centroide, celda_idx in posibles:
                if celda_idx not in usados:
                    centroides.append({
                        "geometry": row.geometry,
                        "centroide": centroide,
                        "zonainfluencia": bloque,
                        "celda_idx": celda_idx
                    })
                    conteo_celdas[celda_idx] += 1
                    usados.add(celda_idx)
                if len(centroides) == 2:
                    break
            # Si no se lograron 2 diferentes, intenta con jitter
            while len(centroides) < 2:
                row, centroide, celda_idx = posibles[len(centroides)]
                # Aplica jitter dentro del polígono si la celda está saturada
                if conteo_celdas[celda_idx] > esperada:
                    for _ in range(20):  # Intenta hasta 20 jitters
                        centroide_jitter = random_point_in_polygon(row.geometry)
                        res_jitter = tree.query(centroide_jitter)
                        if len(res_jitter) > 0:
                            celda_jitter = res_jitter[0]
                            if conteo_celdas[celda_jitter] < conteo_celdas[celda_idx]:
                                centroide = centroide_jitter
                                celda_idx = celda_jitter
                                break
                centroides.append({
                    "geometry": row.geometry,
                    "centroide": centroide,
                    "zonainfluencia": bloque,
                    "celda_idx": celda_idx
                })
                conteo_celdas[celda_idx] += 1
            bloques_centroides.extend(centroides)
        # Calcula el chi2 de esta simulación
        indices = [c["celda_idx"] for c in bloques_centroides]
        conteo = collections.Counter(indices)
        chi2_stat = sum(((conteo.get(i, 0) - esperada) ** 2) / esperada for i in range(num_celdas))
        if chi2_stat < best_chi2:
            best_chi2 = chi2_stat
            best_centroides = bloques_centroides

    # Usa la mejor distribución encontrada
    bloques_centroides = best_centroides
    gdf_grouped_proj = gpd.GeoDataFrame(bloques_centroides, crs=gdf_inter_proj.crs)

    # === OPTIMIZACIÓN LOCAL PARA MEJORAR LA DISTRIBUCIÓN ===
    max_iter = 1000
    improved = True
    iter_count = 0
    best_chi2 = float('inf')
    best_centroides_opt = deepcopy(bloques_centroides)
    num_celdas = len(gdf_grilla_alimentadores_proj)
    num_centroides = len(bloques_centroides)
    esperada = num_centroides / num_celdas

    def calcular_chi2(centroides):
        indices = [c["celda_idx"] for c in centroides]
        conteo = collections.Counter(indices)
        return sum(((conteo.get(i, 0) - esperada) ** 2) / esperada for i in range(num_celdas))

    centroides_actual = deepcopy(bloques_centroides)
    chi2_actual = calcular_chi2(centroides_actual)
    best_chi2 = chi2_actual

    while improved and iter_count < max_iter:
        improved = False
        conteo = collections.Counter([c["celda_idx"] for c in centroides_actual])
        sobre = [i for i in range(num_celdas) if conteo.get(i,0) > esperada]
        bajo = [i for i in range(num_celdas) if conteo.get(i,0) < esperada]
        for idx_s in sobre:
            for c_idx, c in enumerate(centroides_actual):
                if c["celda_idx"] == idx_s:
                    # Intentar mover el centroide dentro de su polígono a una celda menos saturada (jitter)
                    for idx_b in bajo:
                        for _ in range(30):
                            # Generar un punto aleatorio dentro del polígono original
                            poligono = c["geometry"]
                            minx, miny, maxx, maxy = poligono.bounds
                            p_jitter = None
                            for _ in range(10):
                                p = Point(random.uniform(minx, maxx), random.uniform(miny, maxy))
                                if poligono.contains(p):
                                    p_jitter = p
                                    break
                            if p_jitter is None:
                                p_jitter = poligono.representative_point()
                            # Verificar si cae en la celda destino
                            if gdf_grilla_alimentadores_proj.geometry[idx_b].contains(p_jitter):
                                centroides_copia = deepcopy(centroides_actual)
                                centroides_copia[c_idx]["centroide"] = p_jitter
                                centroides_copia[c_idx]["celda_idx"] = idx_b
                                chi2_nuevo = calcular_chi2(centroides_copia)
                                if chi2_nuevo < chi2_actual:
                                    centroides_actual = centroides_copia
                                    chi2_actual = chi2_nuevo
                                    improved = True
                                    break
                        if improved:
                            break
                    if improved:
                        break
            if improved:
                break
        iter_count += 1
    # Si mejoró, actualiza la mejor distribución
    if chi2_actual < best_chi2:
        best_chi2 = chi2_actual
        best_centroides_opt = centroides_actual
    # Usa la distribución optimizada
    bloques_centroides = best_centroides_opt
    gdf_grouped_proj = gpd.GeoDataFrame(bloques_centroides, crs=gdf_inter_proj.crs)

    # Definir la función de reproyección
    project = pyproj.Transformer.from_crs(32717, 4326, always_xy=True).transform
    gdf_grouped_proj["centroide"] = gdf_grouped_proj["centroide"].apply(lambda p: transform(project, p))

    # ====    ========================================================

    # Parroquias
    gdf_rurales = gpd.read_file(GJSON_RURAL).rename(columns={"DPA_DESPAR": "nombre"})
    gdf_urbanas = gpd.read_file(GJSON_URB).rename(columns={"dpa_despar": "nombre"})
    gdf_buses = gpd.read_file(GJSON_BUSES).to_crs("EPSG:4326")
    gdf_metro = gpd.read_file(GJSON_METRO).to_crs("EPSG:4326")

    gdf_rurales["tipo"] = "rural"
    gdf_urbanas["tipo"] = "urbana"
    gdf_parroquias = pd.concat(
        [
            gdf_rurales[["nombre", "geometry", "tipo"]],
            gdf_urbanas[["nombre", "geometry", "tipo"]],
        ],
        ignore_index=True,
    ).set_crs("EPSG:4326")


    # ================== TEST DE CUADRANTES ==================
    indices = [c["celda_idx"] for c in bloques_centroides]
    
    conteo = collections.Counter(indices)
    num_celdas = len(gdf_grilla_alimentadores_proj)
    num_centroides = len(bloques_centroides)

    # Verificaciones de integridad
    print(f"=== VERIFICACIONES DE INTEGRIDAD ===")
    print(f"Total de centroides: {num_centroides}")
    print(f"Centroides contados: {sum(conteo.values())}")
    print(f"Celdas con centroides: {len(conteo)}")
    print(f"Celdas vacías: {num_celdas - len(conteo)}")
    
    # Verificar que todos los centroides fueron contados
    if sum(conteo.values()) != num_centroides:
        print("⚠️  ADVERTENCIA: El número de centroides contados no coincide con el total")
    
    # Frecuencia esperada por celda si la distribución fuera uniforme
    esperada = num_centroides / num_celdas

    # Cálculo del estadístico Chi²
    chi2_stat = sum(((conteo.get(i, 0) - esperada) ** 2) / esperada for i in range(num_celdas))

    # Resultados del test
    df = num_celdas - 1
    valor_critico = chi2.ppf(0.95, df)
    print("=== RESULTADO DEL QUADRANT TEST (2 centroides por alimentador) ===")
    print(f"Número de alimentadores (centroides considerados): {num_centroides}")
    print(f"Número de celdas: {num_celdas}")
    print(f"Frecuencia esperada por celda: {esperada:.2f}")
    print(f"Grados de libertad (df): {df}")
    print(f"Valor crítico Chi² (α=0.05): {valor_critico:.2f}")
    print(f"Estadístico Chi²: {chi2_stat:.2f}")
    if chi2_stat <= valor_critico:
        print("✅ No se rechaza la hipótesis nula (distribución compatible con aleatoriedad)")
    else:
        print("❌ Se rechaza la hipótesis nula (distribución no compatible con aleatoriedad)")

    # ============================================================
    # 🔹 2. CÁLCULO DE GRILLA PARA PARROQUIAS
    # ============================================================

    gdf_parroquias_m = gdf_parroquias.to_crs(epsg=32717)
    gdf_parroquias_m["area_m2"] = gdf_parroquias_m.geometry.area
    mediana_area_parr = gdf_parroquias_m["area_m2"].median()
    lado_celda_parr = mediana_area_parr**0.5
    minx_p, miny_p, maxx_p, maxy_p = gdf_parroquias_m.total_bounds
    grid_cells_p = []
    x = minx_p
    while x < maxx_p:
        y = miny_p
        while y < maxy_p:
            grid_cells_p.append(box(x, y, x + lado_celda_parr, y + lado_celda_parr))
            y += lado_celda_parr
        x += lado_celda_parr

    gdf_grilla = gpd.GeoDataFrame(geometry=grid_cells_p, crs="EPSG:32717").to_crs(
        "EPSG:4326"
    )

    # Intersección celda × parroquia
    gdf_grid_parr = gpd.overlay(
        gdf_grilla[["geometry"]],
        gdf_parroquias[["geometry", "nombre"]],
        how="intersection",
    ).explode(ignore_index=True)

    gdf_grid_parr_proj = gdf_grid_parr.to_crs("EPSG:32717")
    gdf_grid_parr["centroide"] = gdf_grid_parr_proj.centroid.to_crs("EPSG:4326")

    # Mapa
    m = folium.Map(location=[-0.20, -78.50], zoom_start=11, tiles="cartodbpositron")
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

    # Capa de grilla para alimentadores
    fg_grilla_alim = folium.FeatureGroup(name="Grilla Alimentadores", show=False).add_to(m)
    for _, row in gdf_grilla_alimentadores.iterrows():
        folium.GeoJson(
            row.geometry.__geo_interface__,
            style_function=lambda _: {
                "fillColor": "none",
                "color": "brown",
                "weight": 0.5,
                "fillOpacity": 0,
            },
        ).add_to(fg_grilla_alim)

    # Capa de centroides para alimentadores
    fg_centroides_alim = folium.FeatureGroup(
        name="Centroides Alimentadores", show=True
    ).add_to(m)
    for _, row in gdf_grouped_proj.iterrows():
        c = row["centroide"]
        bloque = row["zonainfluencia"]
        folium.CircleMarker(
            location=[c.y, c.x],
            radius=3,
            color="brown",
            fill=True,
            fill_opacity=0.9,
            tooltip = f"Centroide del {bloque}",
        ).add_to(fg_centroides_alim)

    # Estaciones buses
    fg_buses = folium.FeatureGroup(name="Estaciones de Buses", show=False).add_to(m)

    for _, row in gdf_buses.iterrows():
        geom = row.geometry
        folium.GeoJson(
            geom.__geo_interface__,
            style_function=lambda _: {
                "fillColor": "orange",
                "color": "orange",
                "weight": 1.5,
                "fillOpacity": 0.4,
            },
            tooltip="Estación de Bus",
        ).add_to(fg_buses)

        # Icono en el centro del polígono
        centroide = geom.centroid
        folium.Marker(
            location=[centroide.y, centroide.x],
            icon=folium.Icon(icon="bus", prefix="fa", color="orange"),
            tooltip="Estación de Bus",
        ).add_to(fg_buses)

    # Estaciones metro
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

        # Icono centrado
        centroide = geom.centroid
        folium.Marker(
            location=[centroide.y, centroide.x],
            icon=folium.Icon(icon="subway", prefix="fa", color="purple"),
            tooltip=nombre_estacion,
        ).add_to(fg_metro)

    # Paradas de Buses (puntos)
    gdf_paradas = gpd.read_file(GJSON_PARADAS_BUSES).to_crs("EPSG:4326")
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

    # ALIMENTADORES
    fg_alimentadores_padre = folium.FeatureGroup(name="Zonas de Alimentadores").add_to(m)

    # Colores únicos por nombre
    nombre_ids = [
        alimentador_nombre_map.get(aid, aid) for aid in gdf_alimentadores["alimentadorid"].unique()
    ]
    color_map = {
        name: f"#{random.randint(0, 0xFFFFFF):06x}" for name in nombre_ids
    }

    # Crear subcapas por nombre
    subcapas_alimentadores = {}

    for aid, group in gdf_alimentadores.groupby("alimentadorid"):
        nombre = alimentador_nombre_map.get(aid, aid)  # Usa ID si no hay nombre
        if nombre not in subcapas_alimentadores:
            subcapas_alimentadores[nombre] = folium.FeatureGroup(name=nombre).add_to(fg_alimentadores_padre)

        color = color_map[nombre]
        for _, row in group.iterrows():
            folium.GeoJson(
                row.geometry.__geo_interface__,
                style_function=lambda _, col=color: {
                    "fillColor": col,
                    "color": col,
                    "weight": 1,
                    "fillOpacity": 0.4,
                },
                tooltip=nombre,
            ).add_to(subcapas_alimentadores[nombre])

    # Universidades
    df_uni = pd.read_excel(EXCEL_PATH, sheet_name=SHEET_UNI).rename(
        columns=lambda c: c.strip()
    )
    df_uni["UNIVERSIDAD"] = df_uni["UNIVERSIDAD"].str.strip()

    grupo_uni_fin = {"PUBLICA": [], "PRIVADA": []}
    for tipo in ["PUBLICA", "PRIVADA"]:
        fg = folium.FeatureGroup(name=f"Universidades {tipo.title()}").add_to(m)
        grupo_uni_fin[tipo] = fg
        for _, row in df_uni[df_uni["FINANCIAMIENTO"].str.upper() == tipo].iterrows():
            uni = row["UNIVERSIDAD"]
            folium.Marker(
                location=[row["LATITUD"], row["LONGITUD"]],
                title=uni,
                tooltip=f"{uni} – {row['CAMPUS']}",
                icon=folium.Icon(
                    color=(
                        "red"
                        if uni.upper() == "UNIVERSIDAD DE LAS AMERICAS"
                        else "blue"
                    ),
                    icon="university",
                    prefix="fa",
                ),
                careers=[],
            ).add_to(fg)

    # Capa de grilla
    fg_grilla = folium.FeatureGroup(name="Grilla", show=False).add_to(m)
    for _, row in gdf_grilla.iterrows():
        folium.GeoJson(
            row.geometry.__geo_interface__,
            style_function=lambda _: {
                "fillColor": "none",
                "color": "#888",
                "weight": 0.5,
                "fillOpacity": 0,
            },
        ).add_to(fg_grilla)

    # Capa de centroides
    fg_centroides_parr = folium.FeatureGroup(name="Centroides por Intersección", show=False).add_to(m)
    for _, row in gdf_grid_parr.iterrows():
        c = row["centroide"]
        folium.CircleMarker(
            location=[c.y, c.x],
            radius=3,
            color="black",
            fill=True,
            fill_opacity=0.9,
            tooltip="Centroide intersección",
        ).add_to(fg_centroides_parr)

    TreeLayerControl(
        overlay_tree=[
            {
                "label": "Parroquias",
                "layer": fg_parroquias,
            },
            {
                "label": "Universidades",
                "select_all_checkbox": True,
                "children": [
                    {"label": "Públicas", "layer": grupo_uni_fin["PUBLICA"]},
                    {"label": "Privadas", "layer": grupo_uni_fin["PRIVADA"]},
                ],
            },
            {
                "label": "Transporte Público",
                "select_all_checkbox": True,
                "children": [
                    {"label": "Estaciones de Buses", "layer": fg_buses},
                    {"label": "Estaciones de Metro", "layer": fg_metro},
                    {"label": "Paradas de Buses", "layer": fg_paradas},
                ],
            },
            {
                "label": "Grillas",
                "select_all_checkbox": True,
                "children": [
                    {"label": "Grilla Parroquias", "layer": fg_grilla},
                    {"label": "Centroides Parroquias", "layer": fg_centroides_parr},
                    {"label": "Grilla Alimentadores", "layer": fg_grilla_alim},
                    {"label": "Centroides Alimentadores", "layer": fg_centroides_alim},
                ],
            },
            {
                "label": "Alimentadores",
                "select_all_checkbox": True,
                "children": [
                    {"label": nombre, "layer": subcapas_alimentadores[nombre]}
                    for nombre in sorted(subcapas_alimentadores)
                ],
            },
        ]
    ).add_to(m)

    return render_template(
        "index.html",
        mapa=m.get_root().render(),
        map_name=m.get_name(),
        periodos=periodos,
        selected_periodo=selected_periodo,
        ruta_activa="alimentadores"  
    )

