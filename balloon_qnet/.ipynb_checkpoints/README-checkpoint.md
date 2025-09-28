# balloon_qnet

Repository containing NetSquid simulator of balloon-based quantum networks including downlink, uplink and horizontal free-space channels. 

[![DOI](https://zenodo.org/badge/860917433.svg)](https://doi.org/10.5281/zenodo.15388481)

![balloon_net_example](spherical_earth.png)

## Article

The article assosciated to this repository can be currently found on arXiv, under the name [Free-space model for balloon-based quantum networks (arXiv:2412.03356)](https://arxiv.org/abs/2412.03356). 

If you are planning to use this code or data in your research, please cite:

```
@misc{karakostaamarantidou2024freespacemodelballoonbasedquantum,
      title={Free-space model for a balloon-based quantum network}, 
      author={Ilektra Karakosta-Amarantidou and Raja Yehia and Matteo Schiavon},
      year={2024},
      eprint={2412.03356},
      archivePrefix={arXiv},
      primaryClass={quant-ph},
      url={https://arxiv.org/abs/2412.03356}, 
}
```

## Installation 

`balloon_qnet` uses the [NetSquid](https://netsquid.org/) Python package. To install and use NetSquid, you need to first create an account.

For the calculation of the atmospheric transmittance over the link, `balloon_qnet` contains a pre-compiled version of the [`lowtran-piccia` package](https://github.com/francescopiccia/lowtran-piccia).

Please take note that both NetSquid and lowtran-piccia only work on Linux and MacOS. For Windows users it is recommended to use a virtual machine like [Windows Subsystem for Linux (WSL)](https://learn.microsoft.com/en-us/windows/wsl/install). 

For WSL, you can either use directly the terminal (we recommend [Windows Terminal](https://apps.microsoft.com/detail/9n0dx20hk701?hl=en-us&gl=ES)) or the [Visual Studio code extension](https://code.visualstudio.com/docs/remote/wsl). For visualizing plots, you may also install the [`tkinter` package](https://docs.python.org/3/library/tkinter.html):
```
apt install python3-tk
```

### Steps
To install `balloon_qnet` from the source, you first have to clone this repository and then we recomment to create a virtual environment with the appropriate version of Python (the current version of the package has been tested and should be usable for **Python versions 3.12.3 and above**). Specifically, you have to run the following:

```
git clone https://github.com/RajaYehia/balloon_qnet
cd balloon_qnet
python -m venv .venv
source .venv/bin/activate
pip install --extra-index-url https://<username>:<password>@pypi.netsquid.org -e .
```
where `username` and `password` are your NetSquid credentials.

## Structure
This is the main structure of this repository:
```
balloon_qnet/
├── .gitignore
├── LICENSE
├── README.md
├── pyproject.toml
├── spherical_earth.png
├── balloon_qnet/
├── data/
├── Plot files/
└── Studies/
```
In particular:
* `balloon_qnet` contains the source code of the package, including the necessary modules and the compiled version of `lowtran-piccia`.
* `Studies` has the scripts used to generate the different scenarios presented in the paper.
* `Data` contains the data included in the paper, which are the outputs of the scripts in `.txt` format. 
* `Plot files` are the scripts used to generate the plots in the paper for each scenario considered.

## Contributors

Ilektra Karakosta-Amarantidou - University of Padova - [ilektra.karakostaamarantidou@unipd.it](mailto:ilektra.karakostaamarantidou@unipd.it) \
Raja Yehia - ICFO - [raja.yehia@icfo.eu](mailto:raja.yehia@icfo.eu) \
Matteo Schiavon - LIP6 (Sorbonne University) - [matteo.schiavon@lip6.fr](mailto:matteo.schiavon@lip6.fr)


## License

This project is licensed under the MIT License. See the `LICENSE` file for details.