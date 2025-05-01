from flask import Flask, render_template, request
import geopandas as gpd
import pandas as pd
import folium
from shapely.geometry import Point
import branca.colormap as cm
import datetime
from folium.plugins.treelayercontrol import TreeLayerControl
import itertools

app = Flask(__name__)

# ──────────────────────────── RUTAS / ARCHIVOS ────────────────────────────
EXCEL_PATH = "data/universidades_colegios.xlsx"
SHEET_UNI = "Universidades"
SHEET_COL = "Colegios"
CSV_EST = "data/UbicacionEstudiantesPeriodo.csv"
GJSON_RURAL = "data/parroquiasRurales.geojson"
GJSON_URB = "data/parroquiasUrbanas.geojson"
CARRERAS_PATH = "data/baseCarreras.xlsx"


# ───────────────────────────────────────────────────────────────────────────
@app.route("/")
def mapa():
    # 1. ---------------- Periodo seleccionado -----------------
    selected_periodo = request.args.get("periodo")

    # 2. ---------------- Datos de estudiantes -----------------
    df_all = pd.read_csv(CSV_EST, sep=";").rename(columns={"Semestre": "periodo"})
    df_all["periodo"] = df_all["periodo"].astype(str)
    df_all["Latitud"] = pd.to_numeric(df_all["Latitud"], errors="coerce")
    df_all["Longitud"] = pd.to_numeric(df_all["Longitud"], errors="coerce")
    df_all.dropna(subset=["Latitud", "Longitud"], inplace=True)

    periodos = sorted(df_all["periodo"].unique())
    if selected_periodo not in periodos:
        selected_periodo = periodos[0]
    df_est = df_all[df_all["periodo"] == selected_periodo]

    # 3. ---------------- Parroquias y conteo ------------------
    gdf_rurales = gpd.read_file(GJSON_RURAL)
    gdf_urbanas = gpd.read_file(GJSON_URB)
    gdf_rurales["tipo"] = "rural"
    gdf_urbanas["tipo"] = "urbana"
    gdf_rurales = gdf_rurales.rename(columns={"DPA_DESPAR": "nombre"})
    gdf_urbanas = gdf_urbanas.rename(columns={"dpa_despar": "nombre"})
    gdf_parroquias = pd.concat(
        [
            gdf_rurales[["nombre", "geometry", "tipo"]],
            gdf_urbanas[["nombre", "geometry", "tipo"]],
        ],
        ignore_index=True,
    ).set_crs("EPSG:4326")

    gdf_est = gpd.GeoDataFrame(
        df_est,
        geometry=[Point(xy) for xy in zip(df_est["Longitud"], df_est["Latitud"])],
        crs="EPSG:4326",
    )
    gdf_join = gpd.sjoin(gdf_est, gdf_parroquias, how="inner", predicate="within")
    conteo = gdf_join.groupby("nombre").size().reset_index(name="n_estudiantes")
    gdf_parroquias = gdf_parroquias.merge(conteo, on="nombre", how="left")
    gdf_parroquias["n_estudiantes"].fillna(0, inplace=True)
    gdf_parroquias["geometry"] = gdf_parroquias["geometry"].simplify(
        0.0005, preserve_topology=True
    )

    # 4. ---------------- Mapa base ----------------------------
    m = folium.Map(location=[-0.20, -78.50], zoom_start=11, tiles="cartodbpositron")

    # 5. ---------------- Parroquias agrupadas -----------------
    bins = (
        gdf_parroquias["n_estudiantes"]
        .quantile([0, 1 / 3, 2 / 3, 1])
        .round(0)
        .astype(int)
        .tolist()
    )
    gradientes = [
        ["#deebf7", "#9ecae1", "#3182bd"],
        ["#e5f5e0", "#a1d99b", "#31a354"],
        ["#fff7bc", "#fec44f", "#d95f0e"],
    ]
    grupos_parroquia = []
    for i in range(3):
        lwr, upr = bins[i], bins[i + 1]
        sub = gdf_parroquias.query("n_estudiantes >= @lwr & n_estudiantes <= @upr")
        if sub.empty:
            continue
        nombre = f"Grupo {i+1}: {lwr}-{upr}"
        fg = folium.FeatureGroup(name=nombre).add_to(m)
        grupos_parroquia.append({"label": nombre, "layer": fg})
        scale = cm.LinearColormap(gradientes[i], vmin=lwr, vmax=upr)
        for _, row in sub.iterrows():
            folium.GeoJson(
                row.geometry,
                style_function=lambda _, r=row, sc=scale: {
                    "fillColor": sc(r.n_estudiantes),
                    "color": "black",
                    "weight": 0.2,
                    "fillOpacity": 0.7,
                },
                tooltip=f"{row.nombre}: {int(row.n_estudiantes)} estudiantes",
            ).add_to(fg)

    # 6. ---------------- Universidades ------------------------
    df_uni = pd.read_excel(EXCEL_PATH, sheet_name=SHEET_UNI).rename(
        columns=lambda c: c.strip()
    )
    df_uni["UNIVERSIDAD"] = df_uni["UNIVERSIDAD"].str.strip()

    # Carrera ↔ universidad (AÑO 2023)
    df_carr = pd.read_excel(CARRERAS_PATH)
    df_carr = df_carr[df_carr["AÑO"] == 2023]
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
                    color="red" if uni.upper() == "UNIVERSIDAD DE LAS AMERICAS" else "blue",
                    icon="university",
                    prefix="fa",
                ),
                careers=uni_to_carr.get(uni, []),
            ).add_to(fg)

    # 7. ---------------- Colegios por tipo --------------------
    df_col = pd.read_excel(EXCEL_PATH, sheet_name=SHEET_COL).rename(
        columns=lambda c: c.strip()
    )
    df_col["LATITUD"] = pd.to_numeric(df_col["LATITUD"], errors="coerce")
    df_col["LONGITUD"] = pd.to_numeric(df_col["LONGITUD"], errors="coerce")
    df_col.dropna(subset=["LATITUD", "LONGITUD"], inplace=True)

    colegios_grupos, color_cycle = [], itertools.cycle(["orange", "cadetblue"])
    for tipo in sorted(df_col["TIPO"].unique()):
        nombre = f"Colegios {tipo}"
        fg = folium.FeatureGroup(name=nombre).add_to(m)
        colegios_grupos.append({"label": nombre, "layer": fg})
        color = next(color_cycle)
        for _, row in df_col[df_col["TIPO"] == tipo].iterrows():
            folium.Marker(
                location=[row["LATITUD"], row["LONGITUD"]],
                tooltip=row["COLEGIO"].strip(),
                icon=folium.Icon(color=color, icon="graduation-cap", prefix="fa"),
            ).add_to(fg)

    # 8. ---------------- Árbol de capas -----------------------
    overlay_tree = [
        {
            "label": "Parroquias",
            "select_all_checkbox": "Todos",
            "children": grupos_parroquia,
        },
        {
            "label": "Universidades",
            "select_all_checkbox": "Todas",
            "children": [
                {"label": "Públicas", "layer": grupo_uni_fin["PUBLICA"]},
                {"label": "Privadas", "layer": grupo_uni_fin["PRIVADA"]},
            ],
        },
        {
            "label": "Colegios",
            "select_all_checkbox": "Todos",
            "children": colegios_grupos,
        },
    ]

    TreeLayerControl(overlay_tree=overlay_tree, collapsed=False).add_to(m)

    # 9. ---------------- NIVEL → Facultad → Carreras ----------------
    facultades_por_nivel = {}

    for _, r in df_carr.iterrows():
        nivel = r["NIVEL"].strip().upper()
        facultad = r["FACULTAD"].strip()
        carrera = r["CARRERA"].strip()

        if facultad.upper() == "NO COMPARABLE":
            continue

        if nivel not in facultades_por_nivel:
            facultades_por_nivel[nivel] = {}

        if facultad not in facultades_por_nivel[nivel]:
            facultades_por_nivel[nivel][facultad] = set()

        facultades_por_nivel[nivel][facultad].add(carrera)

    # Convertir sets a listas ordenadas
    for nivel in facultades_por_nivel:
        for fac in facultades_por_nivel[nivel]:
            facultades_por_nivel[nivel][fac] = sorted(facultades_por_nivel[nivel][fac])

    # 10. --------------- Render ------------------------------
    return render_template(
        "index.html",
        mapa=m.get_root().render(),  # ⬅⬅ map HTML sin iframe
        map_name=m.get_name(),  # nombre "map_XXXX"
        periodos=periodos,
        selected_periodo=selected_periodo,
        now=datetime.datetime.now(),
        facultades=facultades_por_nivel,
    )


if __name__ == "__main__":
    app.run(debug=True)
