from app import app

app.config['BOOTSTRAP_SERVE_LOCAL'] = True

if __name__=='__main__':
    app.run(debug=False)
