# MDUtils

Utilities for MobiDetails. It is essentially python scripts for MobiDetails administrators. However, one script is designed for the end-user (create_vars_batch.py). See the [wiki](https://github.com/beboche/MDUtils/wiki) for details.


## Installation


requires python 3

### In regular environment

```bash

git clone https://github.com/beboche/MDUtils.git

```

then

```bash

cd MDUtils

```

You need to install some packages:

#### For a full operation

```bash

python3 -m pip install -r requirements.txt

```

#### if you just require the batch creation script:

```bash

python3 -m pip install -r requirements_create_only.txt

```

Done!

### In a virtual environment


```bash

git clone https://github.com/beboche/MDUtils.git
cd MDUtils
python3 -m venv venv
. venv/bin/activate

```

#### For a full operation

```bash

python3 -m pip install -r requirements.txt

```

#### if you just require the batch creation script:

```bash

python3 -m pip install -r requirements_create_only.txt

```


Done!
