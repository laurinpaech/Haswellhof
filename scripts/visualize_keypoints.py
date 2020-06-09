import cv2
import os
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.patches as patches
import numpy as np
import random

# TODO: add command line args 

img_fp = "../images/pumpkin1.jpg" # which image to use
dest = "../plots/" # where to save the visualization
show = True # show plots otherwise just save them
max_kp = None # limit number keypoints shown in plot (None == no limit)
executable = '../build/src/surf'

plot_OpenCV = True
plot_OpenSURF = False # to test against openSurf first generate descriptors with openSurf 
kp_OpenSURF_fp1 = "../images/pumpkin1_OpenSURF_desc.txt"


def draw_kp(img, kps, laplacian=None, show=True, save_fp=None, n_max=None, title=None, rescale=1):
    ax = plt.gca()
    ax.cla()

    if n_max is not None:
        random.seed(4)
        chosen = list(range(0, len(kps)))
        random.shuffle(chosen)
        kps = [kps[i] for i in chosen[:n_max]]
    
    xs = []
    ys = []
    for i, kp in enumerate(kps):
        x, y = kp.pt
        xs.append(x)
        ys.append(y)

        c = 'r' if laplacian is not None and laplacian[i] else 'b'
        circle = plt.Circle(kp.pt, kp.size*rescale, color=c, fill=False, lw=0.5)
        ax.add_artist(circle)

        # rect = patches.Rectangle((kp.pt[0]-kp.size, kp.pt[1]-kp.size), 2*kp.size, 2*kp.size, color='b', fill=False, lw=0.5)
        # ax.add_patch(rect)

    ax.imshow(img1[:,:,[2,1,0]])
    plt.scatter(xs, ys,s=2)
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
    img_fp = os.path.abspath(img_fp)
    kp_fp1 = os.path.abspath(img_fp+'.desc')

    # read images
    img1 = cv2.imread(img_fp) 
    img1_g = cv2.cvtColor(img1, cv2.COLOR_BGR2GRAY)

    # run our implementation
    os.system(f"{executable} {img_fp} {kp_fp1}")
    
    # load txt files written from our implementation
    our1 = np.loadtxt(kp_fp1, dtype=np.float32)
    kp1_our, desc1_our = our1[:,:4], our1[:,4:]

    draw_kp(img1,
            [cv2.KeyPoint(x,y,size) for x,y,size,_ in kp1_our],
            laplacian=[lap for _,_,_,lap in kp1_our],
            n_max=max_kp,
            title="Keypoints from our implementation",
            save_fp=os.path.join(dest,"kps_ours.png"),
            rescale=2.5,
            show=show
           )
    
    if plot_OpenSURF:
        # # run OpenSURF implementation
        # os.system(f"{executable} {img_fp} {kp_fp1}")
        # os.system(f"{executable} {img2_fp} {kp_fp2}")

        # load txt files written from our implementation
        data1 = np.loadtxt(kp_OpenSURF_fp1, dtype=np.float32)
        kp1_OpenSURF, desc1_OpenSURF = data1[:,:4], data1[:,4:]

        draw_kp(img1,
                [cv2.KeyPoint(x,y,size) for x,y,size,_ in kp1_OpenSURF],
                laplacian=[lap for _,_,_,lap in kp1_OpenSURF],
                n_max=max_kp,
                title="Keypoints from OpenSURF implementation",
                save_fp=os.path.join(dest,"kps_OpenSURF.png"),
                rescale=2.5,
                show=show
               )

    if plot_OpenCV:
        # SURF from OpenCV
        surf = cv2.xfeatures2d.SURF_create(nOctaves=3,
                                        hessianThreshold=1150, 
                                        nOctaveLayers=4,
                                        upright=True,
                                        extended=False)
        kp1_openCV, desc1_openCV = surf.detectAndCompute(img1_g,None)
        draw_kp(img1,
                kp1_openCV,
                n_max=max_kp,
                title="Keypoints from OpenCV implementation",
                save_fp=os.path.join(dest,"kps_openCV.png"),
                rescale=0.25,
                show=show
            )

    