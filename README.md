# WHIRLED PEAS

Code for MRI spiral k-space trajectory and matching gradient waveforms.  Gives exact analytic expressions for gradients, k-space, and sampling density compensation.

![Fig 1](https://user-images.githubusercontent.com/116987713/200666042-4684e3e0-f816-44a7-897a-75bff395a1a6.png)

## Code Description
This code is written as a GPI node in python.  To use it,

1.  Go to [gpilab.com](https://gpilab.com/) and become an awesome GPI user

2.  code of interest will begin in the section that starts with:
    ```python
    def compute(self)
    ```
    and you can use nearly all the python code as is.  You will just have to convert the calls that get and pass data,e.g. mostly any line containing the phrase ```self.xxx```

