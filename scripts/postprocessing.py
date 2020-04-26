import cv2
import os
import numpy as np
from matplotlib import pyplot as plt
import numpy as np

# TODO: add command line args 
# TODO: call the c-application given its path and the image paths

max_matches = 100 # limit number matches shown in plot 

img1_fp = "images/img1.ppm"
img2_fp = "images/img2.ppm"
kp_fp1 = "images/kp1_dummy.txt"
kp_fp2 = "images/kp2_dummy.txt"
dest = "./images/"

def visualize_matches(img1, img2, keypoints_1, keypoints_2, matches, 
                        save_fp=None, max_matches=None, show=True):
    height = max(img1.shape[0], img2.shape[0])
    width = img1.shape[1] + img2.shape[1]
    output = np.zeros((height, width, 3), dtype=np.uint8)
    output[0:img1.shape[0], 0:img2.shape[1]] = img1
    output[0:img2.shape[0], img1.shape[1]:] = img2
    if max_matches is not None:
        matches = matches[:max_matches]
        
    for m in matches:
        left = keypoints_1[m.queryIdx][:2].astype(int)
        right = (keypoints_2[m.trainIdx][:2]
                + np.asarray([img1.shape[1],0])).astype(int)
        cv2.line(output, tuple(left), tuple(right), 
                    (0, 255, 0), lineType= cv2.LINE_AA, thickness=2)
    
    plt.figure(figsize=(11,13))
    plt.imshow(output[:,:,[2,1,0]])
    plt.axis('off')
    if save_fp is not None:
        plt.savefig(save_fp)
    if show:
        plt.show()
    else:
        plt.clf()
        
if __name__ == "__main__":
    img1 = cv2.imread(img1_fp) 
    img1_g = cv2.cvtColor(img1, cv2.COLOR_BGR2GRAY)

    # read images
    img2 = cv2.imread(img2_fp) 
    img2_g = cv2.cvtColor(img2, cv2.COLOR_BGR2GRAY)

    # SURF from OpenCV
    sift = cv2.xfeatures2d.SURF_create()
    keypoints_1, descriptors_1 = sift.detectAndCompute(img1_g,None)
    keypoints_2, descriptors_2 = sift.detectAndCompute(img2_g, None)

    keypoints_1 = np.asarray([kp.pt for kp in keypoints_1])
    keypoints_2 = np.asarray([kp.pt for kp in keypoints_2])

    #feature matching
    bf = cv2.BFMatcher(cv2.NORM_L2, crossCheck=True, )
    matches = bf.match(descriptors_1,descriptors_2)

    visualize_matches(img1, img2, keypoints_1, keypoints_2, matches, 
        max_matches=max_matches, 
        save_fp=os.path.join(dest,"matches_openCV.png"))
    
    # create dummy csv to simulate our output (for now)
    csv1 = np.concatenate((keypoints_1, descriptors_1), axis=1)
    csv2 = np.concatenate((keypoints_2, descriptors_2), axis=1)

    np.savetxt(kp_fp1, csv1)
    np.savetxt(kp_fp2, csv2)
    
    # load txt file from our implementation
    keypoints_1, descriptors_1 = np.loadtxt(kp_fp1, dtype=np.float32)[:,:2], np.loadtxt(kp_fp1, dtype=np.float32)[:,2:]
    keypoints_2, descriptors_2 = np.loadtxt(kp_fp2, dtype=np.float32)[:,:2], np.loadtxt(kp_fp2, dtype=np.float32)[:,2:]

    # feature matching
    bf = cv2.BFMatcher(cv2.NORM_L2, crossCheck=True, )
    matches = bf.match(descriptors_1,descriptors_2)

    visualize_matches(img1, img2, keypoints_1, keypoints_2, matches, 
        max_matches=max_matches, 
        save_fp=os.path.join(dest,"matches_ours.png"))