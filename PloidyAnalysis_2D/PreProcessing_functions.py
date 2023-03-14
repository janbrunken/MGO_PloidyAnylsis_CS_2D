import numpy as np
from skimage import io
from natsort import os_sorted
import glob
import os
from tqdm import tqdm

def match_labels(nuc_label_img, mark_label_img):
    """
    Given a nucelar label image and a corresponding marker label image, adjusts the nucelar labels to match the corresponding marker labels.
    """
    mark_label_img=mark_label_img.astype("uint32")
    nuc_label_img=nuc_label_img.astype("uint32")
    overlapping_labels = nuc_label_img * mark_label_img
    matched_nuc_labels = np.divide(overlapping_labels, nuc_label_img, out=np.zeros(overlapping_labels.shape), where=nuc_label_img!=0)
    matched_nuc_labels = matched_nuc_labels.astype("uint16")

    return matched_nuc_labels

def batch_match_labels(nuc_label_dir, mark_label_dir):
    """
    Given the corresponding nuclear and marker label directories, performs batch processing of "match_labels".
    """
    nuc_labels = os_sorted(glob.glob(nuc_label_dir + "*.tif"))
    mark_labels = os_sorted(glob.glob(mark_label_dir + "*.tif"))
    for i in tqdm(range(len(nuc_labels))):
        nuc_label_img = io.imread(nuc_labels[i])
        mark_label_img = io.imread(mark_labels[i])
        out_dir = os.path.join(nuc_labels[i].rsplit("/",2)[0], "labels_matched")
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        matched_nuc_labels = match_labels(nuc_label_img, mark_label_img)
        io.imsave(out_dir + f"/Label_image_matched_Series_{i+1}.tif", matched_nuc_labels, check_contrast=False)

def remove_unmatched_mark(nuc_matched_label_img, mark_label_img):
    """
    Given a matched nucelar label image and a corresponding marker label image, removes excessive marker labels that do not have a matching nuclear label.
    """
    nuc_ul = np.unique(nuc_matched_label_img)
    mark_ul = np.unique(mark_label_img)
    for ul in mark_ul:
        if not ul in nuc_ul:
            mark_mask = (mark_label_img == ul)
            mark_label_img[mark_mask]=0

    return mark_label_img


def batch_remove_unmatched_mark(nuc_matched_label_dir, mark_label_dir):
    """
    Given the corresponding matched nuclear and mark label directories, performs batch processing of "remove_unmatched_mark".
    """
    nuc_matched_labels = os_sorted(glob.glob(nuc_matched_label_dir + "*.tif"))
    mark_labels = os_sorted(glob.glob(mark_label_dir + "*.tif"))
    for n, img in enumerate(tqdm(nuc_matched_labels)):
        nuc_matched_label_img = io.imread(img)
        mark_label_img = io.imread(mark_labels[n])
        out_dir = os.path.join(mark_labels[n].rsplit("/",2)[0], "labels_matched")
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        matched_mark_labels = remove_unmatched_mark(nuc_matched_label_img, mark_label_img)
        io.imsave(out_dir + f"/Label_image_matched_Series_{n+1}.tif", matched_mark_labels, check_contrast=False)

def remove_unmatched_nuc(nuc_matched_label_img, mark_matched_label_img):
    """
    Given a matched nucelar label image and a corresponding marker label image, removes excessive nuclear labels that do not have a matching marker label.
    """
    nuc_ul = np.unique(nuc_matched_label_img)
    mark_ul = np.unique(mark_matched_label_img)
    for ul in nuc_ul:
        if not ul in mark_ul:
            nuc_mask = (nuc_matched_label_img == ul)
            nuc_matched_label_img[nuc_mask]=0
    return nuc_matched_label_img

def batch_remove_unmatched_nuc(nuc_matched_label_dir, mark_matched_label_dir):
    """
    Given the corresponding matched nuclear and marker label directories, performs batch processing of "remove_unmatched_nuc".
    """
    nuc_matched_labels = os_sorted(glob.glob(nuc_matched_label_dir + "*.tif"))
    mark_matched_labels = os_sorted(glob.glob(mark_matched_label_dir + "*.tif"))
    for n, img in enumerate(tqdm(nuc_matched_labels)):
        nuc_matched_label_img = io.imread(img)
        mark_matched_label_img = io.imread(mark_matched_labels[n])
        out_dir = os.path.join(nuc_matched_labels[n].rsplit("/",2)[0], "labels_matched")
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
        matched_nuc_labels = remove_unmatched_nuc(nuc_matched_label_img, mark_matched_label_img)
        io.imsave(out_dir + f"/Label_image_matched_Series_{n+1}.tif", matched_nuc_labels, check_contrast=False)

def reindex_labels(nuc_matched_label_img, mark_matched_label_img):
    """
    Given a matched nucelar label image and a corresponding marker label image, resets the index of both images.
    """
    nuc_ul = np.unique(nuc_matched_label_img)
    template = np.arange(0,len(nuc_ul), 1)
    for n, ul in enumerate(nuc_ul):
        nuc_mask = (nuc_matched_label_img == ul)
        nuc_matched_label_img[nuc_mask]=template[n]

    mark_ul = np.unique(mark_matched_label_img)
    template = np.arange(0,len(mark_ul), 1)
    for n, ul in enumerate(mark_ul):
        mark_mask = (mark_matched_label_img == ul)
        mark_matched_label_img[mark_mask]=template[n]

    return nuc_matched_label_img, mark_matched_label_img

def batch_reindex_labels(nuc_matched_label_dir, mark_matched_label_directories):
    """
    Given the corresponding matched nuclear and marker label directories, performs batch processing of "reindex_labels".
    """
    nuc_matched_labels = os_sorted(glob.glob(nuc_matched_label_dir + "*.tif"))
    mark_matched_labels = os_sorted(glob.glob(mark_matched_label_directories + "*.tif"))
    for n, img in enumerate(tqdm(nuc_matched_labels)):
        nuc_matched_label_img = io.imread(img)
        mark_matched_label_img = io.imread(mark_matched_labels[n])
        nuc_out_dir = os.path.join(nuc_matched_labels[n].rsplit("/",2)[0], "labels_matched_reindexed")
        mark_out_dir = os.path.join(mark_matched_labels[n].rsplit("/",2)[0], "labels_matched_reindexed")
        if not os.path.exists(nuc_out_dir):
            os.mkdir(nuc_out_dir)
        if not os.path.exists(mark_out_dir):
            os.mkdir(mark_out_dir)
        reind_nuc_labels = reindex_labels(nuc_matched_label_img, mark_matched_label_img)[0]
        reind_mark_labels = reindex_labels(nuc_matched_label_img, mark_matched_label_img)[1]
        io.imsave(nuc_out_dir + f"/Label_image_matched_reindexed_Series_{n+1}.tif", reind_nuc_labels, check_contrast=False)
        io.imsave(mark_out_dir + f"/Label_image_matched_reindexed_Series_{n+1}.tif", reind_mark_labels, check_contrast=False)
