import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.tree import DecisionTreeClassifier
from sklearn import metrics
from sklearn.metrics import confusion_matrix
from matplotlib.colors import ListedColormap

def q1(fileIN):
    dat= pd.read_csv(fileIN)
    surved = dat[dat['Survived/Not'] == 1]
    male =   dat[dat['gender'] == 1]

    print ("There are totally " + str(dat.shape[0]) + " entries")
    print ("There are totally " + str(surved.shape[0] ) + " survived.")
    print ("There are totally " + str(male.shape[0] ) + " male passengers.")
    print ("There are totally " + str(dat.shape[0] - male.shape[0] ) + " female passengers.")

def q2(fileIN):
    dat= pd.read_csv(fileIN)
    surved = dat[dat['Survived'] == 1]
    femaleSurv = surved[surved['Sex'] == 'female']
    male =   dat[dat['Sex'] == 'male']

    print ("There are totally " + str(dat.shape[0]) + " entries")
    print ("There are totally " + str(surved.shape[0] ) + " survived.")
    print ("There are totally " + str(male.shape[0] ) + " male passengers.")
    print ("There are totally " + str(dat.shape[0] - male.shape[0] ) + " female passengers.")
    print ("There are totally " + str(femaleSurv.shape[0] ) + " female survived passengers.")

def main():
    data = pd.read_csv('/Users/jyli/myGit/java4lucy/data/input.txt')

 #   print(data.head())
   # decisionTree(data)

##  test female survival rate
    dFemale = data.copy()
    dFemale = dFemale[dFemale['gender'] == 0]

    dFemale = dFemale[dFemale['Age']<9]
    #print(d2.head())

    dFemale = dFemale[dFemale['Sibling']<2.9]
 #   print (dFemale)
    surv = dFemale['Survived/Not']
 #   print(type(surv))
    asn = surv.to_numpy()
 #   print(asn.sum()/asn.size)
    print("Female survival rate  is: " + str(asn.sum()/asn.size))

    #print(data['gender']==0 & data['Age'] < 9)
   # print(data['Age'] < 9)
   # print(data['gender']==0 )
    d2 = data.copy()



    #d2 = d2[d2['gender']==0]
    d2 = d2[d2['gender']==1]
    #print(d2.head())

    d2 = d2[d2['Age']<9]
    #print(d2.head())

    d2 = d2[d2['Sibling']<2.9]
    #print(d2.head())

   # print(d2['Survived/Not'])
    surv = d2['Survived/Not']
    #print(type(surv))

    asn = surv.to_numpy()
    #print (asn)
    print("proportion is: " + str(asn.sum()/asn.size))


    #data <- data[data['gender']==0]
   # data[~data['gender']==0 & data['Age'] < 9]
   # print(data['gender'].head())


def decisionTree (dataIN, dataTest):
    dt = pd.read_csv(dataIN)
    feature_cols = ['gender',  'Age' , 'SiblingAge']
    X = dt.iloc[:,[0,1,2]].values
    y = dt.iloc[:,3].values

    #from sklearn.model_selection import train_test_split
    X_train, X_test, y_train, y_test =  train_test_split(X,y,test_size = 0.25, random_state= 0)

    #print(X_test[0:5,])
    #feature scaling

    sc_X = StandardScaler()
    X_train = sc_X.fit_transform(X_train)
    X_test = sc_X.transform(X_test)

    #from sklearn.tree import DecisionTreeClassifier

    classifier = DecisionTreeClassifier()
    classifier = classifier.fit(X_train,y_train)

    #prediction
    #Accuracy

    y_pred = classifier.predict(X_test)

    #from sklearn import metrics

    print('Accuracy Score:', metrics.accuracy_score(y_test,y_pred))



    ##  test on the know data
    dt2 = pd.read_csv(dataTest)
    feature_cols = ['gender',  'Age' , 'SiblingAge']
    X2 = dt2.iloc[:,[0,1,2]].values
    y2 = dt2.iloc[:,3].values
    #print(X2[0:5,])

    X2 = sc_X.fit_transform(X2)
    #print(X2[0:5,])
    y2_pred = classifier.predict(X2)
    print('Test Accuracy Score:', metrics.accuracy_score(y2,y2_pred))

    ##  Train on full data set
    X_all = sc_X.fit_transform(X)
    classifier2 = DecisionTreeClassifier()
    classifier2 = classifier2.fit(X_all,y)
    y2_all_pred = classifier2.predict(X2)
    print('Test Accuracy Score with all data:', metrics.accuracy_score(y2,y2_all_pred))


if __name__ == '__main__':
  #  print("This is our homework problem\n")
  #  q1('/Users/jyli/myGit/java4lucy/data/input.txt')
    print("\n\nThis is the original training data from kaggle\n")
    q2('/Users/jyli/myGit/java4lucy/data/train.csv')
  #  main()
   # decisionTree ("/Users/jyli/myGit/java4lucy/data/input.txt", "/Users/jyli/myGit/java4lucy/data/output.txt")

