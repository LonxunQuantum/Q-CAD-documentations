PWkit
=====

Download
--------
You can download PWkit from `LongXun github account <https://github.com/LonxunQuantum/>`_.


1. Clone `PWkit` from github 

.. code-block::

   $ git clone git@github.com:LonxunQuantum/pwkit.git

2. Download Python package `pflow`

.. code-block::
   
   $ 

3. Download a conda environment(`pwkit_env.tar.gz`) for PWkit

.. code-block::

   $ 


Installation on MStation
------------------------

We will display how to install `PWkit` on your own PC.

1. Python 的conda环境配置

.. code-block::
   
   # 1. 解压 pwkit_env 的conda环境
   $ mkdir -p pwkit_env
   $ tar -xzf pwkit_env.tar.gz -C pwkit_env

   # 2. 激活 pwkit_env 环境
   $ source pwkit_env/bin/activate

   # 3. 在 pwkit_env 环境下，安装 pflow 库
   (pwkit_env)$ cd pflow; pip install.


2. PWkit 环境变量设置

.. code-block::
   
   # 1. 设置环境变量
   (pwkit_env)$ export PWKIT_ROOT=<Your_pwkit_folder_path>
   (pwkit_env)$ PATH=${PWKIT_ROOT}/bin:$PATH

   # 2. 设置Conda环境 -- 修改 CONDA_PATH
   (pwkit_env)$ vim pwkit.cfg


3. 卸载

.. code-block::

   # 1. 删除 PWkit 文件夹
   $ rm -rf ${PWKIT_ROOT}

   # 2. 删除 PWkit 的配置文件
   $ rm -rf $HOME/.local/pwkit



How to use PWkit on MCloud
--------------------------

.. code-block::

   # 1. 加载环境
   $ module load pwkit/1.0
   # 2. 运行 PWkit
   $ pwkit
   _ ____      ___ __ ___   __ _| |_
   | '_ \ \ /\ / / '_ ` _ \ / _` | __|  website: http://www.lonxun.com
   | |_) \ V  V /| | | | | | (_| | |_   v1.0.0
   | .__/ \_/\_/ |_| |_| |_|\__,_|\__|  pwmat kit Usage: pwkit -h
   |_|