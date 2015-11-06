
# Training and Running the final model
* Make sure that the All_Features.csv file is in the same directory as your code
* Run the R file 'Final_Model.R' 

# Reproduce Report Figures and results (Python and R needed)
#####  1. Create a repository
* Bash "git clone https://github.com/kuancheng/STAT215A-LAB4.git"

#####  2. Install ipython notebook
* Install Anaconda Scientific Package if you don't have  iPython notebook

#####  3. Exploratory data Analysis
* go to "code" folder in your cloned repository
* copy embryo_data.Rda and embryo_imgs.Rda into code folder
* copy raw_fiji folder into code folder
* run lab4_eda.R file 
* read in embryo_data.Rda 
* read in embryo_imgs.Rda
* Now you reproduce the R EDA part
* To extract fiji features to a unified CSV (7.5 GB is needed)
* run proccessing_data.R
* run Random Forest Features Proccessing.ipynb
* a csv All_Features.csv will be generated (~ 5 MBs)
* run 'Lab_4.ipynb' to reproduce the feature selection, model tunning(SVM, Random Forest), and model selection (ROC curves for RF and KNN)
* run 'Lab_4_2.R' to reproduce KNN and weighted KNN models tunning

