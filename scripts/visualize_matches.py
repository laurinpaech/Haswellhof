import cv2
import os
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.patches as patches
import numpy as np
import random

# TODO: add command line args 

img1_fp = "../images/pumpkin1.jpg"
img2_fp = "../images/pumpkin2.jpg"
dest = "../plots/"
show = True # show plots otherwise just save them
max_matches = 100 # limit number matches shown in plot (None == no limit)
executable = '../build/src/surf'

plot_OpenCV = True
plot_OpenSURF = False # to test against openSurf first generate descriptors with openSurf 
kp_OpenSURF_fp1 = "../images/pumpkin1_OpenSURF_desc.txt"
kp_OpenSURF_fp2 = "../images/pumpkin2_OpenSURF_desc.txt"

def visualize_matches(img1, img2, keypoints_1, keypoints_2, matches,
                        save_fp=None, n_max=None, show=True, title=None):
    height = max(img1.shape[0], img2.shape[0])
    width = img1.shape[1] + img2.shape[1]
    output = np.zeros((height, width, 3), dtype=np.uint8)
    output[0:img1.shape[0], 0:img2.shape[1]] = img1
    output[0:img2.shape[0], img1.shape[1]:] = img2
    
    if n_max is not None:
        random.seed(4)
        chosen = list(range(0, len(matches)))
        random.shuffle(chosen)
        matches = [matches[i] for i in chosen[:n_max]]
        
    for m in matches:
        left = keypoints_1[m.queryIdx][:2].astype(int)
        right = (keypoints_2[m.trainIdx][:2]
                + np.asarray([img1.shape[1],0])).astype(int)
        cv2.line(output, tuple(left), tuple(right), 
                    (0, 255, 0), lineType= cv2.LINE_AA, thickness=2)
    
    plt.figure(figsize=(11,13))
    plt.imshow(output[:,:,[2,1,0]])
    plt.axis('off')
    if title is not None:
        plt.title(title)
    if save_fp is not None:
        plt.savefig(save_fp, bbox_inches='tight')
    if show:
        plt.show()
    else:
        plt.clf()
        

if __name__ == "__main__":
    if dest is not None:
        os.makedirs(dest, exist_ok=True)
    # make paths absolute for executable call
    executable = os.path.abspath(executable)
    img1_fp = os.path.abspath(img1_fp)
    img2_fp = os.path.abspath(img2_fp)
    kp_fp1 = os.path.abspath(img1_fp+'.desc')
    kp_fp2 = os.path.abspath(img2_fp+'.desc')

    # read images
    img1 = cv2.imread(img1_fp) 
    img1_g = cv2.cvtColor(img1, cv2.COLOR_BGR2GRAY)
    img2 = cv2.imread(img2_fp) 
    img2_g = cv2.cvtColor(img2, cv2.COLOR_BGR2GRAY)

    # run our implementation
    os.system(f"{executable} {img1_fp} {kp_fp1}")
    os.system(f"{executable} {img2_fp} {kp_fp2}")
    
    # load txt files written from our implementation
    our1 = np.loadtxt(kp_fp1, dtype=np.float32)
    kp1_our, desc1_our = our1[:,:4], our1[:,4:]
    our2 = np.loadtxt(kp_fp2, dtype=np.float32)
    kp2_our, desc2_our = our2[:,:4], our2[:,4:]
    
    # feature matching
    bf = cv2.BFMatcher(cv2.NORM_L2, crossCheck=True)
    matches = bf.match(desc1_our,desc2_our)
    
    visualize_matches(img1, img2, kp1_our[:,:2], kp2_our[:,:2], matches, 
        n_max=max_matches, 
        save_fp=os.path.join(dest,"matches_ours.png"), 
        title="Matched features from our implementation",
        show=show)
    
    if plot_OpenSURF:
        # # run OpenSURF implementation
        # os.system(f"{executable} {img_fp} {kp_fp1}")
        # os.system(f"{executable} {img2_fp} {kp_fp2}")

        # load txt files written from our implementation
        data1 = np.loadtxt(kp_OpenSURF_fp1, dtype=np.float32)
        kp1_OpenSURF, desc1_OpenSURF = data1[:,:4], data1[:,4:]
        data2 = np.loadtxt(kp_OpenSURF_fp2, dtype=np.float32)
        kp2_OpenSURF, desc2_OpenSURF = data2[:,:4], data2[:,4:]

        # feature matching
        bf = cv2.BFMatcher(cv2.NORM_L2, crossCheck=True)
        matches = bf.match(desc1_OpenSURF,desc2_OpenSURF)

        visualize_matches(img1, img2, kp1_OpenSURF[:,:2], kp2_OpenSURF[:,:2], matches, 
            n_max=max_matches, 
            save_fp=os.path.join(dest,"matches_OpenSURF.png"), 
            title="Matched features from OpenSURF implementation",
            show=show)
    
    if plot_OpenCV:
        # SURF from OpenCV
        surf = cv2.xfeatures2d.SURF_create(nOctaves=3,
                                        hessianThreshold=1150, 
                                        nOctaveLayers=4,
                                        upright=True,
                                        extended=False)
        kp1_openCV, desc1_openCV = surf.detectAndCompute(img1_g,None)
        kp2_openCV, desc2_openCV = surf.detectAndCompute(img2_g, None)
        
        kp1_openCV = np.asarray([kp.pt for kp in kp1_openCV])
        kp2_openCV = np.asarray([kp.pt for kp in kp2_openCV])
        
        #feature matching
        bf = cv2.BFMatcher(cv2.NORM_L2, crossCheck=True, )
        matches = bf.match(desc1_openCV,desc2_openCV)
        visualize_matches(img1, img2, kp1_openCV, kp2_openCV, matches, 
            n_max=max_matches, 
            save_fp=os.path.join(dest,"matches_openCV.png"),
            title="Matched features from OpenCV implementation",
            show=show)
        