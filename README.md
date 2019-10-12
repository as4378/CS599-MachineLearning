# CS599-MachineLearning

This repository for sharing the project work done for this course.
In Project1 we implement the least angle regression algorithm for computing the complete lasso path for given input data.
In particular we are computing the lasso path for prostate cancer data in package faraway as shown below:

![Image](https://github.com/as4378/CS599-MachineLearning/blob/master/Project1/Fig2.PNG)

## Citation: 

Hastie et al. Elementes of Statistical Learning, Figure 3.8

The file Project1/lasso.R contains the code that implments the least angle regression algorithm and plots the lasso path.
The make file contains all the necessay commands to run the code and plot the figure. To redo the analysis change to the directory where lasso.R is and just run make on the command line.

This will produce the figure in fig1.pdf. The details are in this file [File](https://github.com/as4378/CS599-MachineLearning/blob/master/Project1/Project1_documentation.pdf)
