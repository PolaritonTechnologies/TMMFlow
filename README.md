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

```

Compile C++ code once more to be on the save side

```console
sudo apt install g++
g++ -shared -o run_filter_stack.so run_filter_stack.cpp -fPIC -fopenmp
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


```

Activate virtual environment

```console
source venv/bin/activate
```

Install python requirements using

```console
pip install -r requirements.txt
```

## Optimization

