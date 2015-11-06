
# Training and Running the Final Model
* Make sure that the FinalData.csv file is in the same directory as your code
* Run the R file 'Final_Model.R' 

# Reproduce Report Figures and results (Python and R needed)
#####  1. Cloning the repository
* Bash "git clone https://github.com/kuancheng/STAT215A-LAB4.git"

#####  2. Install iPython notebook
* Install Anaconda Scientific Package if you don't have  iPython notebook (one easy way to do it)

#####  3. Running the code
* Go to "code" folder in your cloned repository
* Copy embryo_data.Rda and embryo_imgs.Rda into code folder
* Copy raw_fiji folder into code folder
* Read in embryo_data.Rda 
* Read in embryo_imgs.Rda
* Run lab4_eda.R file to reproduce the EDA figures
* Run proccessing_data.R to extract fiji features to a unified CSV (7.5 GB is needed)
* Run Random_Forest_Features_Proccessing.ipynb to generate the super pixel Fiji features as a csv 'All_Features.csv' will be generated (~ 5 MBs)
* Run 'Lab_4.ipynb' to reproduce features selection, model tunning(SVM, Random Forest), and model selection (ROC curves for RF and KNN)
* Run 'Lab_4_2.R' to reproduce KNN and weighted KNN models tunning, need to read in Final_data.csv

