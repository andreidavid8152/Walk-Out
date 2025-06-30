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

mapa_poblacion_parroquias_bp = Blueprint("mapa_calor_poblacion_parroquias", __name__)

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, "..", "data")

# Rutas de archivos
EXCEL_PATH = os.path.join(DATA_DIR, "universidades_colegios.xlsx")
CSV_EST = os.path.join(DATA_DIR, "ubicacionEstudiantesPeriodo.csv")
GJSON_RURAL = os.path.join(DATA_DIR, "parroquiasRurales.geojson")
GJSON_URB = os.path.join(DATA_DIR, "parroquiasUrbanas.geojson")
GJSON_BUSES = os.path.join(DATA_DIR, "estacionesBuses.geojson")
GJSON_METRO = os.path.join(DATA_DIR, "estacionesMetro.geojson")
GJSON_PARADAS = os.path.join(DATA_DIR, "paradasBuses.geojson")


# =========================================================
# 2. RUTA PRINCIPAL DEL MAPA
# =========================================================
@mapa_poblacion_parroquias_bp.route("/mapacalor/poblacion-parroquias")
def mapa():
 return render_template(
        "mapa_calor_poblacionParroquias.html",
        title="Mapa de Calor - Población por Parroquias",
        ruta_activa="poblacion"
)
