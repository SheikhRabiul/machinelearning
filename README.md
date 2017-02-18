# machinelearning Algorithm ( decision-tree learning, C++ implementation)
Implementation of a Machine Learning Algorithm [ decision-tree learning, the  algorithm is taken from the book Artificial Intelligence, A Modern Approach --Russell and Norvig]

This program requires three input file as command line argument as follows:
1. Attribute file: This is the file that contains all attribute name and their possible values. For clarification see attributes.txt file
2. Training file: This is the file based on which the alrogithm develop it's knowledge base and apply that knowledge or learning in classifying data from the provided test file.
3. Test file: This file contains the data that we want to classify.

This was one of my project that I did for my Anomaly and Intrusion Detection System course to detect network intrusion. I used a different dataset for 
this program. For privacy issue, I am not putting that dataset here, instead I am putting a small and easy to comprehend dataset here. 


Windows: 
Compile: g++ dtreelearning.cpp -o dtreelearning.out

Run: dtreelearning.out attributes.txt trainingdata.txt testdata.txt


Linux:
Compile: g++ dtreelearning.cpp -o dtreelearning.out

Run: ./dtreelearning.out attributes.txt trainingdata.txt testdata.txt


#All input data and executable should be in same directory
#The program will also write the output in the same directory as following file name:

output-accuracy-on-training-data.txt

output-accuracy-on-unseen-test-data.txt
