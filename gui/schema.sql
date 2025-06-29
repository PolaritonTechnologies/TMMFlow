DROP TABLE IF EXISTS users;
DROP TABLE IF EXISTS jobs;

CREATE TABLE users (
    id INTEGER PRIMARY KEY,
    username VARCHAR(80) UNIQUE NOT NULL,
    password VARCHAR(128)
);


-- Insert default users
INSERT INTO users (username, password) VALUES ('julian', 'OugaPw');
INSERT INTO users (username, password) VALUES ('florian', 'OugaPw');

CREATE TABLE jobs (
    opt_id INTEGER PRIMARY KEY AUTOINCREMENT,
    job_id INTEGER,
    time_stamp DATETIME DEFAULT CURRENT_TIMESTAMP,
    username VARCHAR NOT NULL,
    filter_name VARCHAR NOT NULL
    optimization_method TEXT,
    current_json TEXT,
    current_data TEXT,
    steps INTEGER,
    initial_merit FLOAT,
    current_merit FLOAT
);