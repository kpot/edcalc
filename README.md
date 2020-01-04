# Electronic Design Calculators

This repository is a collection of jupyter notebooks helpful in building
various DIY electronic devices that require quite a bit of math to make them
right.

Having the design laid out as a jupyter notebook with a bunch of tweakable
parameters can greatly speed up design iterations and give you a refresher later
on the theory of how it all works. As a bonus, by swapping the parameters,
a new device can be easily calculated any time you need.

There are plenty of advanced commercially available tools with similar
and better functionality, but they are not open source and the math behind
them is obscure.

Currently available notebooks:

* [Regulated Half-Bridge Power Supply](./off_line_half_bridge.ipynb)
* [Push-Pull Converter](./off_line_half_bridge.ipynb)
* [Identifier of salvaged magnetic cores](./core_param_estimation.ipynb)


## How to use the notebooks

First, download the repository. On Mac/Linux, assuming you have git installed,
you can do this by running command

```
git clone https://github.com/kpot/edcalc
```
Alternatively, just [download the repository](https://github.com/kpot/edcalc/archive/master.zip) as a zip file and unpack it somewhere.

Next, make sure you have [Python installed](https://wiki.python.org/moin/BeginnersGuide/Download) (version 3.8 or newer).
Then install the libraries you need by launching these commands (from the directory
you have cloned/unpacked the repository into):

```
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

once everything is installed, from the same directory launch

```
.venv/bin/jupyter lab
```

then just open the notebook you need.