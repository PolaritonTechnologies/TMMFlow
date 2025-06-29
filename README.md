# TMM

Simulation and Optimisation webapp based on the Transfer Matrix Method

## Visual Features

#### Interactive Filter Design

<img width="756" alt="Design" src="https://github.com/user-attachments/assets/70017392-2f13-43f1-b48c-ac7eb6c915ba" />

#### Filter Simulation

<img width="927" alt="Simulation" src="https://github.com/user-attachments/assets/eb61cd12-eeaa-439a-b111-8f57028b97d4" />

### Filter Optimisation

<img width="776" alt="Optimisation" src="https://github.com/user-attachments/assets/33cefe9c-c3d8-4b9a-951e-6a40f82c3695" />

## Development Setup

Get a clean ubuntu 22 machine and install gcc and g++ with relevant packages (>
version 11) and git via

```console
sudo apt update
sudo apt install git
sudo apt install
sudo apt install software-properties-common
sudo add-apt-repository ppa:ubuntu-toolchain-r/test
sudo apt update
sudo apt install g++
sudo apt-get install libeigen3-dev
sudo apt-get install nlohmann-json3-dev
```

Compile C++ code

```console
g++ -I/usr/include/eigen3 -I/usr/include/nlohmann -I/usr/local/include -I./include -shared -o interface.so interface.cpp FilterStack.cpp core.cpp input.cpp -fPIC -fopenmp
```

Install python3.10

```console
sudo apt update
sudo apt install software-properties-common
sudo add-apt-repository ppa:deadsnakes/ppa
sudo apt update
sudo apt install python3.10
sudo apt install python3.10-venv
```

Create virtual environment using venv called venv

```console
python3 -m venv venv
```

Activate virtual environment

```console
source venv/bin/activate
```

Install python requirements using

```console
pip install -r requirements.txt
```

## Database

### Initialize DB

Execute from TMM top folder

```console
export PYTHONPATH=$PYTHONPATH:./src
flask --app gui init-db
```

### Add the first user - Configure username and password inside add_user_to_db.py
```console
python3 add_user_to_db.py
```


## Task Scheduling

We're using a redis server that is started on systems with systemctl
automatically after installation (sudo apt-get install redis-server).

On WSL, however, we have to manually start it using "redis-server" and confirm
that it is running using redis-cli ping which should result in the response
"PONG".

In addition to redis server we are using the rq queuing system that is more
lightweight than Celery and easier to use for task scheduling, queuing etc.

Start a default rq worker to actually execute the task within the "src" folder while
having venv activated:

```console
rq worker
```

To start a worker for a specific team, we need to execute

```console
rq worker "<team-name>"
```

instead. The simulations will first look for free workers in the team, and then
default back to default workers. If no workers are available, the job will be queued.

### Job Management

Rq jobs can be managed via the console or in the webbased interface: https://python-rq.org/docs/monitoring/

It would be great if we had implement this monitor in a separate tab at least for admins.

### Structure

database.db

users table

| id  | user    | pw  | e-mail          | team                        | active |
| --- | ------- | --- | --------------- | --------------------------- | ------ |
| 1  | florian | xx  | xx@koeln.de | University of Cologne, HCNB | True   |

jobs table

| opt_id | job_id | time_stamp                | username | filter_name    | optimization_method | initial_json | current_data | steps | initial_merit | current_merit | description                     |
| ------ | ------ | ------------------------- | -------- | -------------- | ------------------- | ------------ | ------------ | ----- | ------------- | ------------- | ------------------------------- |
| 1      | 1      | 2024-07-09 09:28:19.12327 | florian  | test           | None                | {dict}       |              |       |               |               |                                 |
| 2      | 1      | 2024-07-09 09:28:19.12327 | florian   | test           | Nelder-Mead         | {dict}       | {dict}       | 13240 | 13000         | 500           | Optimization with [Nelder-Mead] |
| 3      | 2      | 2024-07-09 09:30:10.12327 | florian   | bandpass_660nm | None                | {dict}       |              |       |               |               |                                 |

- The username is the 1:n connector between the user table and the jobs table.

materials table

| id  | name   | creation_time             | username | team          | material_class | data   |
| --- | ------ | ------------------------- | -------- | ------------- | -------------- | ------ |
| 1   | Ag     | 2024-07-09 09:28:19.12327 | florian   | HCNB, Cologne | default        | {dict} |

### GUI - Main Interaction

The following is how the GUI interacts with the

## Optimization

```

```
