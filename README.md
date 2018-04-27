-----------------------------------------
#          gradient_descent_exp
-----------------------------------------

## Description
 This code allows to fit the data over 
 the regular grid (corresponding to the bcc latice for physical 
 applications) using the adapted gradient descent algorithm
 with the elementary fitting function on each center of the form:
                       f_i(x)=A_i*exp(-B_i*x)+C_i, 
 where A,B and C are the fitting parameters of the center "i". 
 Then, the final fitting function has the form of the sum of elementary functions over the grid:
                        F(x)=sum_i{f_i(x)}

 For the moment the code works for 5 centers with given coordinates 
 which are (ratom(i,:)) but it can be easily generalized to any number 
 and any geometrical configuration of parameters. The code allows precise
 and efficient fitting parameters optimization for different particular
 cases.The code also performs the integration of the fitted function
 over the selected region of 2D space.(however, this part is still 
 in the developement phase).

## Installation 
 The installation needs only the FORTRAN compiler. The default compiler is set to `gfortran`. You can change the compiler in the `Makefile`. To compile the code, just do the `make` command in the `./src` directory. 

## Deinstallation
 To deistall program you can just remove all the file from your computer.
 `make clean` will remove all the compiled files _(which could be useful if you want to be sure that you compile the code from scratch)_

## Usage

`./gradient_descent_exp data_file.dat`

### Technical details
 For the input you need to specify the file name which contains two columns: `x` and `y`. No additional formatting of file is needed.
 The program will output the fitting, integration and linear regression values. Also it will create following files: `fitted.dat`, `param.dat` and `integral.dat`.
- The `fitted.dat` file is the fitted curve.
- The `param.dat` is the file which shows the fitting parameters minimization.
- The `integral.dat` is the integral parameters file.

### Exemple
 You can run this program with the test file `exemple.dat` which you can find in the same directory. 

## Author

Ivan Maliyov

PhD Student
- Doctoral School: EDPIF
- University: Université Paris-Sud, France

