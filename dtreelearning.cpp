/**
* Purpose: Implementation of Decision Tree Learning algorithm
* Author: Sheikh Rabiul Islam , Email: sislam42@students.tntech.edu
* Date: 10/11/2016
*/
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cmath>

using namespace std;

/** Clas to handles file operations */
class IO
{
	public:
	void readFile(std::string &input_file_name,vector <vector <string> > &);	
};

/** Struct to store the tree */
struct treeNode
{
	string divideUpon; /* always column name on which different branch intruduced based on different values of column.
	for leaf node there is no more dividation so it holds value of parents divideUpon.*/
	bool isLeaf; 
	string tag; //Always holds values of a column, usually value of parent, sometimes class level if leaf node. 
	vector <string> childList; // holds list of childs, all distinct values of a column
	vector <treeNode*> child; // holds pointer to each child node
};

/** Helper class for making the tree and testing with test data*/
class Tree
{
	public:
	treeNode* makingDecisionTree_debug(vector <vector <string> > &,vector <vector <string> > &, treeNode*);
	treeNode* makingDecisionTree(vector <vector <string> > &,vector <vector <string> > &,vector <vector <string> > &, treeNode*);
	vector <double> entropy(vector <vector <string> > &,vector <vector <string> > &,vector <vector <string> > &,int, vector <int> &);
	vector <int> count(vector <vector <string> > &,vector <vector <string> > &, vector <vector <string> > &,int);
	vector <int> count2(vector <vector <string> > &,vector <vector <string> > &, vector <vector <string> > &,int,string &);
	void infoGain();
	string importance(vector <vector <string> > &,vector <vector <string> > &,vector <vector <string> > &);	
	void  printDTree(treeNode* node, int);
	bool is_same_class (vector <vector <string> > &);
	bool is_attributes_empty (vector <vector <string> > &,vector <vector <string> > &);
	vector <vector <string> > select_subset_of_data(vector <vector <string> > &, int, string &);
	vector <vector <string> > select_subset_of_attributes(vector <vector <string> > &, string &);
	int col_index_in_attributes_table(vector <vector <string> > &, string &);
	int returnChildSerial(vector <string> &, string &);
	string feedTestDataInToTree(treeNode* nodePtr,vector <string> &line, vector <vector <string> > &attributes, string d_class);
	void accuracy_calculation_and_print(vector <string> &,vector <string> &, string &);
	string majorityValue(vector <vector <string> > &,vector <vector <string> > &);
};


/** Find the best column to split upon based on information gain */
string Tree:: importance(vector <vector <string> > &training_data,vector <vector <string> > &attributes,vector <vector <string> > &attributes_subset)
{	
	vector <double> entropies;
	vector <double> gains;
	vector <int> counts;
	
	int att_index=attributes_subset.size()-1;
	counts=this->count(training_data,attributes,attributes_subset,att_index);
	entropies=this->entropy(training_data,attributes,attributes_subset,att_index,counts);

	for (int i=0;i<training_data[0].size()-1;i++) // for each column except last one 
	{
		//calculate information gain	
		// for each value of the column
		double gain=0.0;
		double entropy_single=0.0;
		
		vector <int> counts3=this->count(training_data,attributes,attributes_subset,i);
		
		for(int j=0;j<attributes_subset[i].size()-1;j++) 
		{	
			vector <int> counts_temp=this->count2(training_data,attributes,attributes_subset,i,attributes_subset[i][j+1]);
			vector <double> entropies_tmp= this->entropy(training_data,attributes,attributes_subset,i,counts_temp);
			double ent_temp = 0.0;
			ent_temp = (double(counts3[j])/double(counts3.back()))*entropies_tmp.back();
			entropy_single +=ent_temp;
		}
		gain = double(entropies.back()-entropy_single);
		gains.push_back(double(entropies.back()-entropy_single));
	}
	

	double max_info_gain=-999999.00;
	int max_info_gain_index=0;
	for(int i=0;i<gains.size();i++)
	{	
		if(gains[i]>max_info_gain)
		{
			max_info_gain=gains[i];
			max_info_gain_index=i;
		}
	}
	
	return attributes_subset[max_info_gain_index][0];
	
}

/** Returns a vector of count, last element of the vector is the sum of count */
vector <int> Tree:: count2(vector <vector <string> > &training_data,vector <vector <string> > &attributes,vector <vector <string> > &attributes_subset,int att_index, string &att_value)
{
	vector <int> counts;
	int sum=0;
	vector <vector <string> > training_data_sub;
	
	for(int ii=0;ii<attributes_subset[attributes_subset.size()-1].size()-1;ii++)
	{
		counts.push_back(0);
	}
	
	for(int i=0;i<training_data.size();i++)
	{	
		if(training_data[i][att_index]==att_value) 
		{
			vector <string> single_line;
			for(int j=0;j<training_data[i].size();j++)
			{ 
					single_line.push_back(training_data[i][j]);
			}	
			if(!single_line.empty()) 
			training_data_sub.push_back(single_line);
		}
		
	}
	
	for(int i=0;i<training_data_sub.size();i++)
	{	
		for(int j=1;j<attributes_subset[attributes_subset.size()-1].size();j++)
		{ 
			if(training_data_sub[i][training_data_sub[i].size()-1]==attributes_subset[attributes_subset.size()-1][j]) 
			{
				counts[j-1]++;
				sum++;
			}
		}		
	}
	counts.push_back(sum);
	return counts; 
}

/** Returns a vector of count, last element of the vector is the sum of count */
vector <int> Tree:: count(vector <vector <string> > &training_data,vector <vector <string> > &attributes,vector <vector <string> > &attributes_subset,int att_index)
{
	vector <int> counts;
	int sum=0;
	for(int ii=0;ii<attributes_subset[att_index].size()-1;ii++)
	{
		counts.push_back(0);
	}
	
	for(int i=0;i<training_data.size();i++)
	{	
		for(int j=1;j<attributes_subset[att_index].size();j++)
		{ 
			if(training_data[i][att_index]==attributes_subset[att_index][j]) 
			{
				counts[j-1]++;
				sum++;
			}
		}		
	}
	counts.push_back(sum);
	return counts; 
}

/** Calculate accuracy, print summary result ont hte screen, write detail result in output file */
void Tree:: accuracy_calculation_and_print(vector <string> &interpolateClass,vector <string> &actual_class, string &title)
{
	int correct_count=0;
	ofstream outFile;
	string file_name="output.txt";
	if(title=="Accuracy of the Tree on test (unseen) data")
	{
		file_name="output-accuracy-on-unseen-test-data.txt";
	}else{
		file_name="output-accuracy-on-training-data.txt";
	}
	
	outFile.open(file_name.c_str());
	// if file opening problem then show message
	if(!outFile)
	{
		cout << "File creation failed for writing output. \n";
	}
	
	outFile << "Detail: " << title << ": "<<endl;
	//outFile << setw(5) << "Row# " << setw(10) << " ActualClass " << setw(10) << "PredictedClass " << setw(10) << " Correct? " << endl<< endl;   
	
	for(int m=0;m<actual_class.size();m++)
	{
		int row=m+1;
		
		if(actual_class[m]==interpolateClass[m])
		{
			correct_count++;
			//outFile << setw(5) << row << setw(10) << actual_class[m] << setw(10) << interpolateClass[m] << setw(10) << " ---Yes--- " << endl; 
			outFile << "Row: " << row << endl;
			outFile << "Actual Class: " << actual_class[m] << endl;
			outFile << "Predicted Class: " << interpolateClass[m] << endl;
			outFile << "Correct?: " << "Yes" << endl;
			outFile << endl;
			
		}else{
			//outFile << setw(5) << row << setw(10) << actual_class[m] << setw(10) << interpolateClass[m] << setw(10) << " ----No---- " << endl; 
			
			outFile << "Row: " << row << endl;
			outFile << "Actual Class: " << actual_class[m] << endl;
			outFile << "Predicted Class: " << interpolateClass[m] << endl;
			outFile << "Correct?: " << "No" << endl;
			outFile << endl;
		}
	}
	
	int total_records=actual_class.size();
	double accuracy = (double(correct_count)/double(total_records))*100;
		
	outFile << endl;
	cout <<endl << title << ":" <<endl;
	cout << "Total records: " << total_records << endl;
	cout << "Correctly predicted class: " << correct_count << endl;
	cout << "Accuracy:  " << accuracy << " % "<< endl;
	cout << "Detail result is written in to file: " <<  file_name << endl;
	
	outFile<< "Brief: " << title << ":" << endl;
	outFile << "Total records: " << total_records << endl;
	outFile << "Correctly predicted class: " << correct_count << endl;
	outFile << "Accuracy:  " << accuracy << " % "<< endl;
	outFile.close();
}

/** Recursively print the decision tree */
void Tree:: printDTree(treeNode* node,int level_counter)
{    
	if(node == NULL) {
		return;
	}		
	string tab_="\t";
	string newline_="\n";
	if (!node->child.empty()) 
	{  
		for(int i=1;i<=level_counter;i++)
		{
			tab_=tab_+ "\t";
			newline_=newline_+ "\n";
		}
		cout << newline_<< "Parent: " << node->tag << tab_ <<" Splitting Column: " << node->divideUpon <<endl;
		
		int iii;
		cout  << endl;
		cout << tab_ <<" childrens:";
		for (iii = 0; iii < node->childList.size(); iii++) {   
			cout << "  " << node->childList[iii];
		}
		cout << endl <<endl;
		for (iii = 0; iii < node->child.size(); iii++) {   
			this->printDTree(node->child[iii],level_counter);
		}
		return;
    } else {
		cout << tab_ <<" class=" << node->tag<<" , ";
		return;
	}
		level_counter++;
		cout <<endl;
		cout <<endl;
}


/** Calculate entropy*/
vector <double> Tree:: entropy(vector <vector <string> > &training_data,vector <vector <string> > &attributes, vector <vector <string> > &attributes_subset, int att_index,vector <int> &counts)
{	
	vector <double> entropies;
	double entropy_tot = 0.00;
	
	for(int i=0;i<attributes_subset[attributes_subset.size()-1].size()-1;i++)
	{
		double entropy_tmp=0.00;
		
		if(counts[i]!=0 && counts.back()!=0)  // as (0/10) log (0/10) will give errror 
		{
			entropy_tmp=-((double(counts[i])/double(counts.back())*(log(double(counts[i])/double(counts.back()))))/log(2));
		}
		entropy_tot+=entropy_tmp;
		entropies.push_back(entropy_tmp);
	}
	entropies.push_back(entropy_tot);
	return entropies;
	
}

/** Calculate most frequent class */
string Tree:: majorityValue(vector <vector <string> > &training_data,vector <vector <string> > &attributes_subset)
{
	vector <int> counts;
	int att_index=attributes_subset.size()-1;
	for(int ii=0;ii<attributes_subset[att_index].size()-1;ii++)
	{
		counts.push_back(0);
	}
	
	for(int i=0;i<training_data.size();i++)
	{	
		for(int j=1;j<attributes_subset[att_index].size();j++)
		{ 
			if(training_data[i][att_index]==attributes_subset[att_index][j]) 
			{
				counts[j-1]++;
			}
		}		
	}
	
	int max_count=-999999;
	int max_count_index=0;
	
	for(int i=0;i<counts.size();i++)
	{	
		if(counts[i]>max_count)
		{
			max_count=counts[i];
			max_count_index=i;
		}
	}
	
	return attributes_subset[att_index][max_count_index+1];
 
}

/** Return serial number of the column in the table */
int Tree:: col_index_in_attributes_table(vector <vector <string> > &attributes, string &col_name)
{	
	for(int t=0;t<attributes.size();t++)
	{
		if(attributes[t][0]==col_name) return t;
	}
}

/** Return child serial number among all child of the node*/
int Tree:: returnChildSerial(vector <string>  &childList,string &child_name)
{	
	for(int t=0;t<childList.size();t++)
	{
		if(childList[t]==child_name) return t;
	}
	
}

/** feed one row of test data to the tree and return predicted class based on the traverse result of the tree */
string Tree :: feedTestDataInToTree( treeNode* node,vector <string> &line, vector <vector <string> > &attributes, string d_class)
{
	string class_predict;
	while (!node->isLeaf && !node->child.empty()) 
	{
		int col_index = col_index_in_attributes_table(attributes,node->divideUpon);
		string value = line[col_index];
		
		int childSerial = returnChildSerial(node->childList, value);
		node = node->child[childSerial];
		
		if (node == NULL)
		{
			class_predict = d_class;
			break;
		}
		
		class_predict = node->tag;
	}
	
	return class_predict;
}

/** Select subset of data applicable for a particular branch */
vector <vector <string> > Tree:: select_subset_of_data(vector <vector <string> > &training_data,int col, string &col_value)
{
	vector <vector <string> > training_data_subset;
	
	for (int i = 0; i <training_data.size(); i++) { 
		
		vector <string> single_row;
		
		if(training_data[i][col]==col_value)
		{	
			for (int j = 0; j < training_data[i].size(); j++) { 
				if(j!=col)
				{
					single_row.push_back(training_data[i][j]);
				}
			} 
		}
		
		if(!single_row.empty())
		training_data_subset.push_back(single_row);
	}
	return training_data_subset;
}

/** Select subset of attrbutes remaining */
vector <vector <string> > Tree:: select_subset_of_attributes(vector <vector <string> > &attributes, string &col_name)
{
	vector <vector <string> > attributes_subset;
	
	for (int i = 0; i < attributes.size(); i++) { 
		vector <string> single_row;
		if(attributes[i][0]!=col_name)
		{	
			for (int j = 0; j < attributes[i].size(); j++) { 
				
				single_row.push_back(attributes[i][j]);
			} 
		}
		if(!single_row.empty())
		attributes_subset.push_back(single_row);
	}		
	return attributes_subset;
}

/** Return true if all remaining rows are from same class */
bool Tree:: is_same_class (vector <vector <string> > &training_data)
{
	for(int i=0;i<training_data.size();i++)
	{
		if(i!=training_data.size()-1)
		{
			if(training_data[i][training_data[0].size()-1]!=training_data[i+1][training_data[0].size()-1]) return false;
		}else{
			if(training_data[i][training_data[0].size()-1]!=training_data[0][training_data[0].size()-1]) return false;
		}
	}
	
	return true;
}

/** Return true if there is still some training data to process but no attribute left (except last one) for processing, this may happen due to noise in data */
bool Tree:: is_attributes_empty (vector <vector <string> > &training_data,vector <vector <string> > &attributes_subset)
{
	if(attributes_subset.size()<=1 && training_data.size()>0 )
	{
		return true;
	}	
	return false;
}

/** Recursively make the decision tree and return the final tree */
treeNode* Tree:: makingDecisionTree(vector <vector <string> > &training_data,vector <vector <string> > &attributes,vector <vector <string> > &attributes_subset, treeNode*  node)
{
	if(training_data.size() == 0)   // no more rows left for processing in training data set
	{
		return NULL;
		
	}else if(this->is_same_class(training_data)){  // all remaining rows in training set are from same class
		
		node->tag=training_data[0][training_data[0].size()-1];
		node->isLeaf=true;
		
	}else if(this->is_attributes_empty(training_data,attributes_subset)){	// still some rows left in training set but attribute to process is empty (except last one). This may happen due to noise in data	
		node->tag=this->majorityValue(training_data, attributes_subset);
		node->isLeaf=true;
	}else{	
		string divideUpon_=importance(training_data,attributes,attributes_subset); // most important attribute to split on based on infromation gain
		node->divideUpon=divideUpon_; 
		int divideUponColIndex=	col_index_in_attributes_table(attributes,divideUpon_); // return column serial number in main attrubute set
		int divideUponColIndexForSubset = col_index_in_attributes_table(attributes_subset,divideUpon_); // return column serial number in  attrubute sub set, we are excluding already processed column
		
		for(int c=1;c<attributes[divideUponColIndex].size();c++)
		{			
			//for each value a new branch
			treeNode* childNode = (treeNode*) new treeNode;
			childNode->divideUpon = divideUpon_;
			childNode->tag = attributes[divideUponColIndex][c];
			childNode->isLeaf = false;
			node->childList.push_back(attributes[divideUponColIndex][c]);
			
			vector <vector <string> > training_data_subset=select_subset_of_data(training_data,divideUponColIndexForSubset,attributes[divideUponColIndex][c]);
			vector <vector <string> > attributes_subset=select_subset_of_attributes(attributes,divideUpon_);

			node->child.push_back(makingDecisionTree(training_data_subset,attributes,attributes_subset,childNode));
		}				
	}
	
	return node;
}

/** Read files provided in commandline and store them in a vecto */
void IO:: readFile(std::string &input_file_name, vector <vector <string> > &input_data)
{
		ifstream inputFile;
		inputFile.open(input_file_name.c_str());
		// if file opening problem then show message
		if(!inputFile)
		{
			cout << "File  opening problem \n";
		}
		
		while (inputFile)
		{
			string single_line;
			// copy whole line in single_line variable
			if (!getline( inputFile, single_line )) break; 
			istringstream sll( single_line );
			vector <string> tuple;

			while (sll)
			{
			  string single_line;
			  // copies untill spaces in to sll variable
			  if (!getline( sll, single_line, ' ' )) break;
			  // insert the column
			  tuple.push_back( single_line );
			}
			// insert the row
			input_data.push_back( tuple );
	  }
		
	  if (!inputFile.eof())
	  {
		cout << "something wrong with the file.\n";
	  }
	  inputFile.close(); // close the file
	  
}

/** Main starts here*/	
int main(int argc, const char *argv[])
{	
	
	// check number of commandline arguments
	if(argc!=4)
	{
		cout << endl << "Please provide three command line arguments and re-run the program. First argument is the attribute file, Second argument is the training file, Third argument is the test file." << endl;
		exit(1);
	}
	
	// creating object of IO class
	IO ioObj;
	
	// reading attribute file into attributes vector of vector
	vector <vector <string> > attributes;
	string attribute_file=argv[1];
	ioObj.readFile(attribute_file, attributes);
	
	//reading training data into training_data vector of vector
	vector <vector <string> > training_data;
	string training_data_file=argv[2];
	ioObj.readFile(training_data_file, training_data);
	
	//reading test data into test_data vector of vector
	vector <vector <string> > test_data;
	string test_data_file=argv[3];
	ioObj.readFile(test_data_file, test_data);
	
	// start making decision tree
	Tree tr; // object of Tree class
	treeNode* rootNode = new treeNode; // instance of structure to hold tree. 
	
	vector <vector <string> > attributes_subset=attributes; // we need to maintain two attribute table, one is original one and another one is excluding one attribute at a time.Initially both are same.  
	rootNode=tr.makingDecisionTree(training_data,attributes,attributes_subset,rootNode); // make the decision tree and return the tree
	
	cout << "Printing Decision Tree: " << endl;
	int level_counter=1;
	tr.printDTree(rootNode,level_counter); // print the decision tree
	cout << endl;
	
	string d_class=tr.majorityValue(training_data,attributes); // find out the most frequeent class and set it as default class.
	
	// start predicting class  by traversing the tree 
	vector <string>  interpolateClass;
	for(int l=0;l<test_data.size();l++) 
	{	// pass one line of test data at a time
		interpolateClass.push_back(tr.feedTestDataInToTree( rootNode,test_data[l], attributes, d_class));
	}
	
	// take actual class in test data 
	vector <string> actual_class;
	for(int kk=0;kk<test_data.size();kk++)
	{
		actual_class.push_back(test_data[kk][test_data[0].size()-1]);
	}
	
	// accuracy on test set by compairing predicted class vs actual class of test data
	string title="Accuracy of the Tree on test (unseen) data";
	tr.accuracy_calculation_and_print(interpolateClass,actual_class,title);
	
	
	// clear container 	
	interpolateClass.clear();	
	actual_class.clear();
	
	// accuracy on training set by compairing predicted class vs actual class of training data
	for(int kk=0;kk<training_data.size();kk++)
	{
		actual_class.push_back(training_data[kk][training_data[0].size()-1]);
	}
	for(int l=0;l<training_data.size();l++)
	{
		interpolateClass.push_back(tr.feedTestDataInToTree( rootNode,training_data[l], attributes, d_class));
	}
    title="Accuracy of the Tree on training data";
	tr.accuracy_calculation_and_print(interpolateClass,actual_class,title);
	
	// Brief result will be printed on the screen and detail result will be written on two file.
	// we are done with the program
	
	return 0;
}