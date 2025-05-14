from flask import Flask
from routes.main import main_bp
from routes.universidades import universidades_bp

app = Flask(__name__)
app.register_blueprint(main_bp)
app.register_blueprint(universidades_bp)

if __name__ == "__main__":
    app.run(debug=True, host="0.0.0.0", port=5000)
