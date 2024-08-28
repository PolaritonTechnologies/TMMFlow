# TMM

Simulation and Optimisation software based on the Transfer Matrix Method

## Development Setup

Get a clean ubuntu 22 machine and install gcc and g++ (> version 11) and git via

```console
sudo apt update
sudo apt install git
sudo apt install
sudo apt install software-properties-common
sudo add-apt-repository ppa:ubuntu-toolchain-r/test
sudo apt update
sudo apt install gcc-11 g++-11
```

%% Not Supported yet %%

Provide correct git credentials and subsequently clone repository (subsequently for an AWS instance)

```
Follow the instructions to create a personal access token if unavailable:
https://docs.github.com/en/get-started/getting-started-with-git/about-remote-repositories#cloning-with-https-urls
and:
[https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/managing-your-personal-access-tokens](https://docs.github.com/en/authentication/keeping-your-account-and-data-secure/managing-your-personal-access-tokens#creating-a-personal-access-token-classic

```

%%%% For a remote instance: %%%%

Install Filezilla on your local computer: https://filezilla-project.org/
Use the key pair of your AWS instance to configure the SFTP connection
Use your instance Public IPv4 address in settings - configure on port 22
Upload the TMM folder to the instance
(In case you need to run quickly, copy only examples,materials,src,tests, requirements)

%%%%

````

Compile C++ code once more to be on the save side

```console
sudo apt install g++
g++ -shared -o run_filter_stack.so run_filter_stack.cpp -fPIC -fopenmp
````

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

### Structure

database.db

users table

| id  | user    | pw  |
| --- | ------- | --- |
| 1   | julian  | xx  |
| 2   | florian | xx  |

jobs table

| opt_id | job_id | time_stamp                | username | filter_name    | optimization_method | initial_json | current_data | steps | initial_merit | current_merit |
| ------ | ------ | ------------------------- | -------- | -------------- | ------------------- | ------------ | ------------ | ----- | ------------- | ------------- |
| 1      | 1      | 2024-07-09 09:28:19.12327 | julian   | test           | None                | {dict}       |              |       |               |               |
| 2      | 1      | 2024-07-09 09:30:10.12327 | julian   | bandpass_660nm | Nelder-Mead         | {dict}       | {dict}       | 13240 | 13000         | 500           |

- The job_id is the 1:n connector between the user table and the jobs table.

## Task Scheduling

We're using a redis server that is started on systems with systemctl
automatically after installation (sudo apt-get install redis-server).

On WSL, however, we have to manually start it using "redis-server" and confirm
that it is running using redis-cli ping which should result in the response
"PONG".

In addition to redis server we are using the rq queuing system that is more
lightweight than Celery and easier to use for task scheduling, queuing etc.

Start an rq worker to actually execute the task within the "src" folder while
having venv activated:

```console
rq worker
```

Also I think it has to be reloaded when the code is changed similarly to the
flask application itself to work.

### Job Management

Rq jobs can be managed via the console or in the webbased interface: https://python-rq.org/docs/monitoring/

It would be great if we had implement this monitor in a separate tab at least for admins.

### GUI - Main Interaction

The following is how the GUI interacts with the

## Optimization

### Design Strategies (2 materials)

1. Start from a generic design (e.g. 25 pairs of SiO2 - Ta2O5 on one side and 5
   layers of MgF2 - Ta2O5 on the back) with thicknesses defined by the approximate wavelength region of interest (lambda/4).
2. Most important point: Define targets to match **all** design features of
   interest. This might also include less important regions but ideally the whole
   spectrum of interest is covered. Consider lower tolerances for critical areas.
3. Run **dual annealing** (~3 h) + 5 x **Nelder-Mead** (~10 min) + **basin
   hopping** (~3 h) to obtain a global best solution.
4. Refine targets to match desired features and repeat steps 1 - 3.
