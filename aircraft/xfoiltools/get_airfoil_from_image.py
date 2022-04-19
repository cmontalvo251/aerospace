import cv2
from scipy.spatial import Delaunay
import numpy as np 
import matplotlib.pyplot as plt

img = cv2.imread('number3.jpg')
#cv2.imshow("original",img)
#cv2.waitKey(0)
print(img.shape)

plt.imshow(img)

x1 = 2500
x2 = 2800
y1 = 1300
y2 = 2000

cropped = img[x1:x2,y1:y2,:]

plt.figure()
plt.imshow(cropped)

edges = cv2.Canny(cropped,50,120)

plt.figure()
plt.imshow(edges)

x12 = 78
x22 = 678
y12 = 130
y22 = 230

edges_cropped = edges[y12:y22,x12:x22]

plt.figure()
plt.imshow(edges_cropped)

plt.show()