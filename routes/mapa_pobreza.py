# ────────────────────────────────────────────────────────────────
#  Flask |  Mapa de Pobreza (NBI) a nivel parroquial + conteo de alimentadores
# ────────────────────────────────────────────────────────────────
from flask import Blueprint, render_template
from folium.plugins.treelayercontrol import TreeLayerControl
import geopandas as gpd
import folium
import pandas as pd
import os, json, unicodedata

mapa_pobreza_bp = Blueprint("mapa_pobreza", __name__)

# ── Rutas de archivos ───────────────────────────────────────────
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, "..", "data")

GJSON_ALIMENTADORES = os.path.join(DATA_DIR, "alimentadores.geojson")
GJSON_RURAL = os.path.join(DATA_DIR, "parroquiasRurales.geojson")
GJSON_URB = os.path.join(DATA_DIR, "parroquiasUrbanas.geojson")
NBI_PATH = os.path.join(DATA_DIR, "NBI.xlsx")  # ajusta si es .csv
ID_ALIM_JSON = os.path.join(DATA_DIR, "idAlimentadores.json")


# ── Normalizar texto (quitar tildes) ────────────────────────────
def normalize(txt: str) -> str:
    txt = unicodedata.normalize("NFKD", str(txt))
    return "".join(c for c in txt if not unicodedata.combining(c)).strip().upper()


# ── Rangos y colores EXACTOS ────────────────────────────────────
RANGE_COLORS = [
    (45.9, 69.2, "#922546"),
    (32.4, 45.8, "#e73b3d"),
    (19.0, 32.3, "#fd9d58"),
    (11.3, 18.9, "#fede8a"),
    (4.8, 11.2, "#ffffd3"),
]


def color_for(val):
    """Devuelve el color correspondiente al % de pobreza por NBI."""
    if pd.isna(val):
        return "#cccccc"
    for low, high, col in RANGE_COLORS:
        if low <= val <= high:  # inclusivo en ambos extremos
            return col
    return "#cccccc"


# ── Diccionario código → nombre de alimentador ──────────────────
with open(ID_ALIM_JSON, encoding="utf-8") as f:
    alimentador_nombre_map = {
        item["code"]: item["name"] for item in json.load(f)["codedValues"]
    }

# ── Colores por bloque de zona de influencia ────────────────────
COLORES_BLOQUES = {
    "Bloque 1": "#a2635c", ##
    "Bloque 2": "#568d9d", ##
    "Bloque 3": "#7d9455", ##
    "Bloque 4": "#775277", ##
    "Bloque 5": "#a68155", ##
    "Bloque 6": "#a39b57",##
    "Bloque 7": "#716f80",
    "Bloque 8": "#a68e9b", ##
    "Bloque 9": "#698d88",##
    "Bloque 10": "#899a87", ##
    "Bloque 12": "#a2a27b",##
    "Bloque 13": "#734e73" ##
}


# ────────────────────────────────────────────────────────────────
@mapa_pobreza_bp.route("/mapa/pobreza")
def mapa():
    # 1. Cargar geojson de capas básicas
    gdf_alim = gpd.read_file(GJSON_ALIMENTADORES).to_crs("EPSG:4326")
    gdf_rural = gpd.read_file(GJSON_RURAL).rename(columns={"DPA_DESPAR": "nombre"})
    gdf_urb = gpd.read_file(GJSON_URB).rename(columns={"dpa_despar": "nombre"})

    gdf_rural["tipo"] = "rural"
    gdf_urb["tipo"] = "urbana"

    # Unión de parroquias en un solo GeoDataFrame
    gdf_parr = pd.concat(
        [
            gdf_rural[["nombre", "geometry", "tipo"]],
            gdf_urb[["nombre", "geometry", "tipo"]],
        ],
        ignore_index=True,
    ).set_crs("EPSG:4326")

    # 2. Leer Excel NBI y unir con parroquias
    df_nbi = pd.read_excel(NBI_PATH, dtype=str)

    df_nbi["NBI_pct"] = (
        df_nbi["Personas con pobreza por NBI (%)"]
        .str.replace("%", "", regex=False)
        .str.replace(",", ".", regex=False)
        .astype(float)
        .round(1)  # 🔑 redondear a 1 decimal
    )
    df_nbi["Parroquia_norm"] = df_nbi["Parroquia"].apply(normalize)

    gdf_parr["nombre_norm"] = gdf_parr["nombre"].apply(normalize)

    # Merge pobreza NBI ↔ parroquias
    gdf_parr = gdf_parr.merge(
        df_nbi[["Parroquia_norm", "NBI_pct"]],
        left_on="nombre_norm",
        right_on="Parroquia_norm",
        how="left",
    )

    # 3. Conteo de alimentadores por parroquia  ──────────────────
    # Spatial join: asignar a cada alimentador la parroquia que intersecta
    alim_parr = gpd.sjoin(
        gdf_alim[["geometry"]],
        gdf_parr[["nombre", "geometry"]],
        how="left",
        predicate="intersects",
    )

    # Serie: parroquia → nº de alimentadores
    conteo_alim = (
        alim_parr.groupby("nombre", observed=True)
        .size()
        .reindex(gdf_parr["nombre"])  # conserva orden original
        .fillna(0)
        .astype(int)
    )

    # 4. Mapa base Folium ────────────────────────────────────────
    m = folium.Map(location=[-0.20, -78.50], zoom_start=11, tiles="OpenStreetMap")
    
    # Agregar CSS para oscurecer el mapa y hacerlo gris medio-oscuro
    custom_css = """
    <style>
    .leaflet-container {
        background-color: #474749 !important;
    }
    .leaflet-tile-pane {
        filter: brightness(0.25) contrast(1.4) saturate(0.1) invert(0.1);
    }
    </style>
    """
    m.get_root().html.add_child(folium.Element(custom_css))

    # 4-a. Capa Alimentadores coloreados por bloque ---------------
    fg_alim = folium.FeatureGroup(name="Alimentadores", show=True).add_to(m)

    for _, row in gdf_alim.iterrows():
        nombre_alim = alimentador_nombre_map.get(
            row["alimentadorid"], row["alimentadorid"]
        )
        
        # Obtener color según zona de influencia (bloque)
        zona_influencia = row.get("zonainfluencia", "")
        color_bloque = COLORES_BLOQUES.get(zona_influencia, "#6f6f70")  # gris por defecto
        
        folium.GeoJson(
            row.geometry.__geo_interface__,
            style_function=lambda _, color=color_bloque: {
                "fillColor": color,
                "color": color,
                "weight": 1.5,
                "fillOpacity": 0.7,
            },
            tooltip=f"{nombre_alim}<br>Zona: {zona_influencia}",
        ).add_to(fg_alim)

    # 4-b. Capa Parroquias coloreadas por NBI --------------------
    fg_parr = folium.FeatureGroup(name="Parroquias – NBI", show=True).add_to(m)

    for _, row in gdf_parr.iterrows():
        pct = row["NBI_pct"]
        pct_txt = f"{pct:.1f} %" if pd.notnull(pct) else "Sin dato"

        folium.GeoJson(
            {
                "type": "Feature",
                "geometry": row.geometry.__geo_interface__,
                "properties": {"Parroquia": row["nombre"], "NBI_pct": pct_txt},
            },
            style_function=lambda _: {
                "fillColor": "transparent",
                "color": "white",
                "weight": 1.0,
                "fillOpacity": 0.0,
            },
            tooltip=folium.GeoJsonTooltip(
                fields=["Parroquia", "NBI_pct"],
                aliases=["Parroquia:", "Pobreza NBI (%):"],
            ),
        ).add_to(fg_parr)

    # 6. Leyenda «Alimentadores por parroquia» (inferior-derecha) --
    legend_rows_alim = [
        f"{parr}: <b>{n}</b>"
        for parr, n in conteo_alim.sort_values(ascending=False).items()
        if n > 0  # omite las de 0 (opcional)
    ]

    legend_alim_html = (
        '<div style="position: fixed; bottom: 30px; right: 30px; z-index:9999;'
        "background:white;padding:8px 10px;border:2px solid gray;border-radius:4px;"
        'font-size:13px;line-height:18px;max-height:280px;overflow-y:auto;">'
        "<b>Alimentadores por Parroquia</b><br>"
        + "<br>".join(legend_rows_alim)
        + "</div>"
    )
    m.get_root().html.add_child(folium.Element(legend_alim_html))

    # 7. Control de capas en árbol --------------------------------
    TreeLayerControl(
        overlay_tree=[
            {"label": "Alimentadores", "layer": fg_alim},
            {"label": "Parroquias – NBI", "layer": fg_parr},
        ]
    ).add_to(m)

    # 8. Render ---------------------------------------------------
    return render_template(
        "mapa_pobreza.html",
        mapa=m.get_root().render(),
        map_name=m.get_name(),
        ruta_activa="pobreza",
    )
