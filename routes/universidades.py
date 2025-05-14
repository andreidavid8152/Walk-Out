from flask import Blueprint, render_template, request
import geopandas as gpd
import pandas as pd
import folium
import datetime
from shapely.geometry import Point
from folium.plugins.treelayercontrol import TreeLayerControl
import os
from utils.helpers import darken_color

universidades_bp = Blueprint("universidades", __name__)

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, "..", "data")

EXCEL_PATH = os.path.join(DATA_DIR, "universidades_colegios.xlsx")
SHEET_UNI = "Universidades"
CSV_EST = os.path.join(DATA_DIR, "ubicacionEstudiantesPeriodo.csv")
GJSON_RURAL = os.path.join(DATA_DIR, "parroquiasRurales.geojson")
GJSON_URB = os.path.join(DATA_DIR, "parroquiasUrbanas.geojson")
CARRERAS_PATH = os.path.join(DATA_DIR, "baseCarreras.xlsx")


@universidades_bp.route("/universidades")
def mapa_universidades():
    # 🎯 Incluye solo lógica de parroquias + universidades + filtros

    # Periodo
    selected_periodo = request.args.get("periodo")
    df_all = pd.read_csv(CSV_EST, sep=";").rename(columns={"Semestre": "periodo"})
    df_all["periodo"] = df_all["periodo"].astype(str)
    periodos = sorted(df_all["periodo"].unique())
    if selected_periodo not in periodos:
        selected_periodo = periodos[0]

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

    # Universidades
    df_uni = pd.read_excel(EXCEL_PATH, sheet_name=SHEET_UNI).rename(
        columns=lambda c: c.strip()
    )
    df_uni["UNIVERSIDAD"] = df_uni["UNIVERSIDAD"].str.strip()
    df_carr = pd.read_excel(CARRERAS_PATH)
    df_carr["PERIODO"] = df_carr["PERIODO"].astype(str)
    if selected_periodo in ["202410", "202420"]:
        df_carr = df_carr[df_carr["PERIODO"].isin([selected_periodo, "202400"])]
    else:
        df_carr = df_carr[df_carr["PERIODO"].isin([selected_periodo, "202520"])]
    uni_to_carr = df_carr.groupby("UNIVERSIDAD")["CARRERA"].apply(list).to_dict()

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
                careers=uni_to_carr.get(uni, []),
            ).add_to(fg)

    TreeLayerControl(
        overlay_tree=[
            {"label": "Parroquias", "layer": fg_parroquias},
            {
                "label": "Universidades",
                "children": [
                    {"label": "Públicas", "layer": grupo_uni_fin["PUBLICA"]},
                    {"label": "Privadas", "layer": grupo_uni_fin["PRIVADA"]},
                ],
            },
        ]
    ).add_to(m)

    facultades_por_nivel = {}
    for _, r in df_carr.iterrows():
        nivel = r["NIVEL"].strip().upper()
        facultad = r["FACULTAD"].strip()
        carrera = r["CARRERA"].strip()
        if facultad.upper() == "SIN REGISTRO":
            continue
        facultades_por_nivel.setdefault(nivel, {}).setdefault(facultad, set()).add(
            carrera
        )
    for nivel in facultades_por_nivel:
        for fac in facultades_por_nivel[nivel]:
            facultades_por_nivel[nivel][fac] = sorted(facultades_por_nivel[nivel][fac])

    return render_template(
        "universidades.html",
        mapa=m.get_root().render(),
        map_name=m.get_name(),
        periodos=periodos,
        selected_periodo=selected_periodo,
        facultades=facultades_por_nivel,
    )
