`Back to manual </docs/manual/README.md>`__

Python package management
=========================

In general, it is best to create a virtual environment before installing
the python packages you need for a certain project. This way, the
package versions required for compatibility with one project don’t
interfere with the compatibility of another project. You will only need
to create the environment and install the packages once. In subsequent
sessions, you just need to activate the environment you previously made.

Setting up a virtual environment
--------------------------------

+------------------+-----------------------------+---------------------+
| Command          | Pip                         | Conda               |
+==================+=============================+=====================+
| Creating an      | ``$ pyth                    | ``$ cond            |
| environment      | on3 -m venv [path/to/env]`` | a create -n [env]`` |
+------------------+-----------------------------+---------------------+
| Activating an    | ``$ source [                | ``$ con             |
| environment      | path/to/env]/bin/activate`` | da activate [env]`` |
+------------------+-----------------------------+---------------------+
| Deactivating an  | ``$ deactivate``            | ``$                 |
| environment      |                             |  conda deactivate`` |
+------------------+-----------------------------+---------------------+

Managing packages
-----------------

+--------------+---------------------------+---------------------------+
| Command      | Pip                       | Conda                     |
+==============+===========================+===========================+
| Installing   | ``$ pip3 install [p       | ``$ conda install [p      |
| packages     | ackage1] [package2] ...`` | ackage1] [package2] ...`` |
+--------------+---------------------------+---------------------------+
| List of      | ``$ pip3 list``           | ``$ conda list``          |
| installed    |                           |                           |
| packages     |                           |                           |
+--------------+---------------------------+---------------------------+

`Other useful terminal commands <./terminal.md>`__
--------------------------------------------------
