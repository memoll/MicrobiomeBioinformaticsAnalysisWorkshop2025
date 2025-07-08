###############################################################
# Microbiome Bioinformatic Analysis                           #
# Introduction to R                                           #
# Author: ArrietaLab - University of Calgary                  #
# Date: July 2025                                             #
# Location: IUCBC - Cordoba (Argentina)                       #
###############################################################

#Let's R                

# The # symbol in R indicates "comments" and is NOT processed by R.
# Why are comments important? They help others - and your future self - understand your code.

#enter to run the code
output <- "Hello World!"
output 

#brackets [] in front of the output show at which value of the output you are situated

#Errors and Warnings ####
#warning: Warns you about a potential issue in the input/output but still executes the function
x <- c("2", -3, "end", 0, 4, 0.2)
as.numeric(x)

#error: Notifies you of a problem that must be fixed to proceed.
x * 10

#Google and ChatGPT are your best friends!

#Object #### 
#see slides
y <- 7
y
A <- "toto"
a <- "foo"
A
a

#Operators #### 
#see slides
#__Arithmetic Operators ####
5+1
10-2
3*5
18/3
2^5
2**5
z <- 14/5
z
y + z

#Now it's your turn
2 + 20 * 31 - 11
#pay attention to the order of operataions
(2 + 20) * (31 - 11)

#built-in constant 
pi 
pi * 2^3

?Constants #in console
month.name

#__Logical Operators ####
y > z
y < z
y + 5 > z

#Functions ####
#see slides

ls
?ls
?median
help(ls)

ls()
mean()
mean(z)

?seq
MySequence <- seq(1, 100, by=5)
MySequence #print your new created sequence

MyMean <- mean(MySequence)
MyMean

SumMySequence <- sum(MySequence)
SumMySequence

print(SumMySequence)


#Let's create some functions
#see slides
my_function <- function() { # create a function with the name my_function
  print("Hello World!")
}
#What is the result of running the line below?
my_function()
#Good jab! You created a function (that does take any arguments)

#Now try this one:
fname_function <- function(fname) {
  paste(fname, "Griffin")
}

fname_function("Peter")
fname_function("Lois")
fname_function("Stewie")
fname_function()
#What does the fname_function do? Does it accept any arguments?
#What happens if you run it empty?

#Let's try another one:
city_function <- function(city = "Tehran") {
  paste("I am from", city)
}

city_function("Shiraz")
city_function("Mashhad")
city_function() # will get the default value, which is Tehran
city_function("Tabriz")


#There are many functions built in to R.
#We can also load sets of functions and data (called packages) to do different types of analysis. 
#We'll talk about packages later.

#R tip: Use the 'Up' and 'Down' arrows in the console to navigate through previous commands.

#Types of data in R ####
?typeof
?class

#__Numeric ####
z
class(z)
typeof(z)

u <- 5.0; typeof(u)
v <- 5; typeof(v)
w <- 5L; typeof(w)

#__Character ####
typeof(toto)
typeof("toto")

c <- 5
d <- "5"
class(c)
class(d)

#__Logical ####
typeof(y > z)
class(y > z)

#__Factor ####
days <- factor(c("Saturday", "Monday", "Monday", "Friday", "Saturday"))
days
class(days)

days[1]
days[1:3]

#Types of objects in R ####
#see slides

#__Vectors ####
weekdays <- c("Monday", "Tuesday", "Wednesday", "Thursday", "Friday")
weekdays 
class(weekdays)
is.vector(weekdays)

weekdays[1]
weekdays[1:3]

names(weekdays) <- c("day1", "day2", "day3", "day4", "day5")
weekdays
weekdays[1:3]
length(weekdays)
dim(weekdays)

#Matrix ####
#numeric matrix
mat <- matrix(c(1, 2, 3, 4, 5, 6), nrow = 2, ncol = 3)
mat

rownames(mat) <- c("row1", "row2")
colnames(mat) <- c("col1", "col2", "col3")
mat
mat[,2]
mat[1,]
mat[2,3]

class(mat)
is.matrix(mat)
dim(mat)

#character matrix
animals <- matrix(c("cat", "dog", "duck", "dog", "cat", "cow"), nrow = 3, ncol = 2)
animals

rownames(animals) <- c("row1","row2","row3")
colnames(animals) <- c("group1","group2")
animals
dim(animals)
animals[3,2]

#Dataframe ####
df <- as.data.frame(mat)
df
df$col2

df$col4 <- c("foo", "toto")
dim(df)
df$col4

class(df)
is.data.frame(df)

# obtain a summary of each column of a data.frame
summary(df)

#List ####
list <- list(foo = weekdays, matrix1 = mat, matrix2 = animals, dataframe = df)
names(list)
list
list$matrix2

class(list)
is.list(list)

#Manage your workspace ####
# what is the current directory?
getwd()

#save workspace image to a file in current directory
save.image("day2_workspace.RData")

#delete all objects in the workspace
rm(list = ls())
ls()

#load the workspace we just saved
load("day2_workspace.RData")
ls()

#Importing data ####
mydata <- read.csv2("PATH_TO_FILE/FILE_NAME.csv", header=TRUE, sep=";", row.names = 1) #use column 1 as row names 
mydata <- read.table("PATH_TO_FILE/FILE_NAME.tsv", header=TRUE, sep="\t") 

#Exporting data ####
write.csv2(mydata, "PATH_TO_FILE/FILE_NAME.csv")

#R tip: Tab autocompletion
#If you type part of an object or function name and hit tab, you'll see a list of objects with names matching that text.
write_

#Data Visualization ####

#load data
data(iris) #built-in data
#explore data
dim(iris)
class(iris)
summary(iris)
head(iris)
colnames(iris)

#histogram: visualize distribution of a dataset
hist(iris$Sepal.Width) #of of sepal width

#boxplot of petal length by species
boxplot(Petal.Length ~ Species, data = iris, xlab = "Species", ylab = "Petal length")

#what type of data (object) it takes?
?boxplot #read the Arguments section

#R Packages ####
#see slides

#Install R packages from one of the following repositories####

#__CRAN:
install.packages("package-name", dependencies = TRUE)

#__Bioconductor:
install.packages("BiocManager")
BiocManager::install("package-name")

#__GitHub:
install.packages("remotes")  
remotes::install_github("author/package-name")
#Or
install.packages("devtools")  
devtools::install_github("author/package-name")

#Load the package once it is installed ####
library(package-name)

#Read the package's documentation ####
help(package-name)

#Install the following R packages for the next session, using one of the above methods (whichever works for you):

#dada2
#phyloseq
#vegan
#ggplot2
#tidyverse 
#dplyr
#BiocGenerics
#SummarizedExperiment
#rstatix
#DESeq2
#ggpubr
#RColorBrewer
#gridExtra

