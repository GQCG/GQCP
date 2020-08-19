---
id: developer_vscode
title: VScode setup
sidebar_label: VScode setup
---

At GQCG our code editor of choice is [Visual Studio code](https://code.visualstudio.com). VScode allows for easy collaboration, uses built-in git commands and offers a lot of useful extensions. In order to use VScode to collaborate on GQCP, some setup is required. This will be covered in this small step-by-step guide. 

## Step 1: Install the needed VScode extensions

To comfortably work on GQCP using VScode, it is recommended you install several extensions.

- __C/C++__: [C/C++ language support](https://marketplace.visualstudio.com/items?itemName=ms-vscode.cpptools)
- __CMake__: [CMake language support](https://marketplace.visualstudio.com/items?itemName=twxs.cmake)
- __CMake Tools__: [Allows for compiling/building within VScode](https://marketplace.visualstudio.com/items?itemName=gocarlos.cmake-tools)


## Step 2: Clone the repository & install the necessary dependencies

How to do this is already covered in [Getting started](getting_started.md). 
You will need the GQCP folder to set up your VScode workspace. By installing the provided `gqcp_dev` environment, the necessary dependencies and their paths will al be set correctly for the following steps.


## Step 3: Setting up your C++ compiler

We recommend using Clang 7.1.0 as your default C++ compiler. On mac, this can easily be installed using [macports](https://www.macports.org). After the installation, you can add the compiler to VScode by editing the `cmake-tools-kits.json` file. (Which can be found by searching for `Edit user-local CMake kits`).

<!--DOCUSAURUS_CODE_TABS-->

<!--json-->
```json
{
    "name": "Clang 7.1.0",
    "compilers": {
      "C": "/opt/local/bin/clang-mp-7.0",
      "CXX": "/opt/local/bin/clang++-mp-7.0"
    }
  }
```

<!--END_DOCUSAURUS_CODE_TABS-->

> **Note:** The paths to your compilers may be different. What's shown here is the most common path when macports is used to install Clang 7.1.0.


After this is done, you can select Clang 7.1.0 as your default compiler in VScode.


## Step 4: Set up your GQCP workspace

After opening the GQCP folder in VScode, select `file, Save Workspace As..`, and save your workspace. This will create a folder `.vscode` with a file `name.code-workspace` in it. Within this file, complete the settings as follows:

<!--DOCUSAURUS_CODE_TABS-->

<!--json-->
```json
{
	"folders": [
		{
			"path": ".."
		}
	],
	"settings": {
		"cmake.configureSettings":{
			"CMAKE_PREFIX_PATH": "${CONDA_PREFIX}/envs/gqcp_dev",
			"CMAKE_CXX_COMPILER":"/opt/local/bin/clang++-mp-7.0",
            "BUILD_TESTS": "TRUE",
            "BUILD_PYTHON_BINDINGS": "TRUE",
            "CMAKE_INSTALL_PREFIX": "${CONDA_PREFIX}/envs/gqcp_dev",
            "PYTHON_EXECUTABLE": "${CONDA_PREFIX}/envs/gqcp_dev/bin/python"
		},
		"files.associations": {
			"random": "cpp",
			"__bit_reference": "cpp",
			"algorithm": "cpp"
		},
		"cmake.cmakePath": "${CONDA_PREFIX}/envs/gqcp_dev/bin/cmake",
		"C_Cpp.default.includePath": [
			"${CONDA_PREFIX}/envs/gqcp_dev/include",
			"${CONDA_PREFIX}/envs/gqcp_dev/include/eigen3",
			"${workspaceFolder}/gqcp/include"
		],
	}
}
```

<!--END_DOCUSAURUS_CODE_TABS-->

## Step 5: Building locally

By using the `CMake tools` extension in your side bar, you should now be able to build and reconfigure GQCP locally. This can be done by clicking on the extension in your sidebar and the clicking on the `...` in the upper right corner. Then you can choose to either `Clean rebuild` or `Clean reconfigure`. 