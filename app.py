from flask import Flask
from routes.main import main_bp
from routes.mapa_calor_universidades import mapa_uni_bp
from routes.mapa_calor_colegios import mapa_colegios_bp
from routes.mapa_calor_empresas import mapa_empresas_bp
from routes.mapa_calor_estudiantes import mapa_estudiantes_bp
from routes.mapa_calor_poblacion_parroquias import mapa_poblacion_parroquias_bp
from routes.mapa_calor_hospitales import mapa_hospitales_bp

app = Flask(__name__)
app.register_blueprint(main_bp)
app.register_blueprint(mapa_uni_bp)
app.register_blueprint(mapa_colegios_bp)
app.register_blueprint(mapa_empresas_bp)
app.register_blueprint(mapa_poblacion_parroquias_bp)
app.register_blueprint(mapa_estudiantes_bp)
app.register_blueprint(mapa_hospitales_bp)

if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=5001)
