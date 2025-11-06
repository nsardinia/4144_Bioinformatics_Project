import numpy as np

from sklearn.neighbors import KNeighborsClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC

def knn(x_train: np.ndarray, x_test: np.ndarray, y_train: np.ndarray) -> np.ndarray:
    model = KNeighborsClassifier()
    y_pred = model.fit(x_train, y_train).predict(x_test)
    return y_pred

def lr(x_train: np.ndarray, x_test: np.ndarray, y_train: np.ndarray) -> np.ndarray:
    model = LogisticRegression()
    y_pred = model.fit(x_train, y_train).predict(x_test)
    return y_pred

def gnb(x_train: np.ndarray, x_test: np.ndarray, y_train: np.ndarray) -> np.ndarray:
    model = GaussianNB()
    y_pred = model.fit(x_train, y_train).predict(x_test)
    return y_pred

def svc(x_train: np.ndarray, x_test: np.ndarray, y_train: np.ndarray) -> np.ndarray:
    model = SVC()
    y_pred = model.fit(x_train, y_train).predict(x_test)
    return y_pred