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

Provide correct git credentials and subsequently clone repository (subsequently for an AWS instance)

```console
# Generate ssh key pair on AWS instance
ssh-keygen -t rsa -b 4096 -C "your_email@example.com"
# Display key
cat ~/.ssh/id_rsa.pub
# Copy displayed public key and add to github account in settings --> SSH and GPG keys --> New SSH key --> paste into key field and press add ssh key
# Clone git repository to the AWS instance
git clone git@github.com:username/repository.git
```

Compile C++ code once more to be on the save side

```console
g++ -shared -o run_filter_stack.so run_filter_stack.cpp -fPIC
```

Install python3.10

```console
sudo apt update
sudo apt install software-properties-common
sudo add-apt-repository ppa:deadsnakes/ppa
sudo apt update
sudo apt install python3.10
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

## Optimization

Optimization benchmark (test_optimisation_DBRCavity.json on Julian's laptop
varying all layer thicknesses plus the organic's position )

| Method                 | Time to Optimize | Comments                                                                         |
| ---------------------- | ---------------- | -------------------------------------------------------------------------------- |
| minimize (Nelder-Mead) | 60 s             | Using 1e-9 tolerance                                                             |
| dual_annealing         | 1802 s           | Using 1e-5 tolerance, 7000 iterations                                            |
| basinhopping           | 126 s            | Does not support bounds (recovers results of minimize for the example)           |
| direct                 | 777 s            | merit only minimized to 0.003 after 2831 calls                                   |
| differential_evolution |                  | no sign of convergence after ~50000 iterations                                   |
| shgo                   |                  | did not even start to optimize                                                   |
| brute                  |                  | not really an option because it requires too many iterations for relevant stacks |

dual_annealing, optimized weights:
[ 43.44974405 32.76934884 31.45234875 190.87887203 65.24474902
140.11114023 19.47791629 132.36720361 324.52413575 81.47214849
69.120849 239.70597544 62.72602507 282.03428307 447.86895182
459.8029801 143.59800317 484.79319159 178.5618158 24.86932182
221.55959766 15.84151904 2.68046251]

![Dual annealing optimized image](/doc/dual-annealing-optimized.png)
