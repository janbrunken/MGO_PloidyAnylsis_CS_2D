import numpy as np
from skimage import io, measure, color
from skimage.measure import perimeter
import matplotlib.pyplot as plt
import pandas as pd
import os
import plotly.graph_objects as go
import plotly.express as px
from tqdm import tqdm


def measure_properties(nuc_image, mark_image, cc_image, nuc_label, mark_label, mark, cc):
    """
    Given a grayscale DNA staining image (e.g. DAPI), its corresponding label image, a grayscale cell type marker image, its corresponding label image
    and a grayscale cell cycle marker image, returns the DNA, cell type marker and cell cycle marker median intensity and integrated density as well as
    the nuclear and marker area for each label and for the background.
    """
    unique_labels = np.unique(nuc_label)
    properties = []
    for ul in unique_labels:
        # Select only pixels belonging to a specific label
        nuc_mask = (nuc_label == ul)
        mark_mask = (mark_label == ul)
        nuc_label_pixels = nuc_image[nuc_mask]
        mark_label_pixels = mark_image[mark_mask]
        cc_label_pixels = cc_image[nuc_mask]

        # Calculate median intensity, integrated density, area and perimeter
        nuc_median_intensity = np.median(nuc_label_pixels)
        nuc_integrated_density = np.sum(nuc_label_pixels)
        nuc_area = np.count_nonzero(nuc_mask)
        nuc_perimeter = perimeter(nuc_mask, neighborhood=0) 

        mark_median_intensity = np.median(mark_label_pixels)
        mark_integrated_density = np.sum(mark_label_pixels)
        mark_area = np.count_nonzero(mark_mask)

        cc_median_intensity = np.median(cc_label_pixels)
        cc_integrated_density = np.sum(cc_label_pixels)

        properties.append({
            'label': ul,
            'nuc_median_intensity': nuc_median_intensity,
            'nuc_integrated_density': nuc_integrated_density,
            'nuc_area': nuc_area,
            'nuc_perimeter': nuc_perimeter,
            f'{mark}_median_intensity': mark_median_intensity,
            f'{mark}_integrated_density': mark_integrated_density,
            f'{mark}_area': mark_area,
            f'{cc}_median_intensity': cc_median_intensity,
            f'{cc}_integrated_density': cc_integrated_density
        })
    return properties

def process_data(data, mark, cc):
    """
    Given results from the measure_properties function, calculates background-corrected cell cycle marker and DNA median intensities,
    cell cycle marker and DNA corrected total cellular fluorescence values and nuclear volume.
    """
    df = pd.DataFrame(data)

    df['nuc_median_intensity_bgCorr'] = df['nuc_median_intensity']-df.iloc[0]['nuc_median_intensity'] #calculate median intensity
    df['nuc_CTCF'] = df['nuc_integrated_density']-df['nuc_area']*df.iloc[0]['nuc_median_intensity'] #calculate CTCF
    df["nuc_volume"] = 4/3*np.pi*np.power(np.sqrt(df["nuc_area"]/np.pi), 3) #calculate nuclear volume
    df['nuc_circularity'] = 4*np.pi*df['nuc_area']/np.power(df['nuc_perimeter'], 2)


    df[f'{mark}_median_intensity_bgCorr'] = df[f'{mark}_median_intensity']-df.iloc[0][f'{mark}_median_intensity'] #calculate median intensity
    df[f'{mark}_CTCF'] = df[f'{mark}_integrated_density']-df[f'{mark}_area']*df.iloc[0][f'{mark}_median_intensity'] #calculate CTCF
    df[f"{mark}_volume"] = 4/3*np.pi*np.power(np.sqrt(df[f"{mark}_area"]/np.pi), 3) #calculate nuclear volume

    df[f'{cc}_median_intensity_bgCorr'] = df[f'{cc}_median_intensity']-df.iloc[0][f'{cc}_median_intensity'] #calculate median intensity
    df[f'{cc}_CTCF'] = df[f'{cc}_integrated_density']-df['nuc_area']*df.iloc[0][f'{cc}_median_intensity'] #calculate CTCF

    df.drop(df.head(1).index,inplace=True) #remove background line
    df['nuc_median_intensity_bgCorr_n'] = df['nuc_median_intensity_bgCorr']/df['nuc_median_intensity_bgCorr'].median() #calculate normalized median intensity
    df['nuc_CTCF_n'] = df['nuc_CTCF']/df['nuc_CTCF'].median() #calculate normalized CTCF

    df[f'{mark}_median_intensity_bgCorr_n'] = df[f'{mark}_median_intensity_bgCorr']/df[f'{mark}_median_intensity_bgCorr'].median() #calculate normalized median intensity
    df[f'{mark}_CTCF_n'] = df[f'{mark}_CTCF']/df[f'{mark}_CTCF'].median() #calculate normalized CTCF

    df[f'{cc}_median_intensity_bgCorr_n'] = df[f'{cc}_median_intensity_bgCorr']/df[f'{cc}_median_intensity_bgCorr'].median() #calculate normalized median intensity
    df[f'{cc}_CTCF_n'] = df[f'{cc}_CTCF']/df[f'{cc}_CTCF'].median() #calculate normalized CTCF

    return df

def process_images(nuc_images, mark_images, cc_images, nuc_label_images, mark_label_images, mark, cc):
    """
    Given input directories for DNA staining images, marker images and cell cycle marker images as well as the corresponding nuclear and marker label images,
    performs measure_properties and process_data for all contained images.
    """

    dfs = pd.DataFrame()
    for i in tqdm(range(len(nuc_images))): # if the last image series is a stitched image, use "for i in range(len(cc_images)-1):"
        nuc_image = io.imread(nuc_images[i])
        mark_image = io.imread(mark_images[i])
        cc_image = io.imread(cc_images[i])
        nuc_labels = io.imread(nuc_label_images[i])
        mark_labels = io.imread(mark_label_images[i])
        mark = mark
        cc = cc
        properties = measure_properties(nuc_image, mark_image, cc_image, nuc_labels, mark_labels, mark, cc)
        df = process_data(properties, mark, cc)
        df["series"] = i+1
        dfs = pd.concat([dfs, df], axis=0, ignore_index = True)
    return dfs

def simple_overlay(image, labels):

    """
    Given an inupt image and its corresponding label image, displays the input image as well as a masked overlay.
    """

    plt.rcParams['figure.dpi'] = 150

    img=io.imread(image)
    mask=io.imread(labels)

    ratio = np.amax(img) / 256
    img = (img/ratio).astype("uint8")

    ratio = np.amax(mask) / 256
    mask = (mask/ratio).astype("uint8")

    fig = plt.figure(figsize=(12,8))
    fig.add_subplot(1, 2, 1)
    plt.imshow(img, cmap="Greys_r")
    plt.axis('off')
    plt.title(f"Image")

    fig.add_subplot(1, 2, 2)
    io.imshow(color.label2rgb(mask, img, bg_label=0))
    plt.axis('off')
    plt.title(f"Masked image")



def interactive_overlay(image_files, labels_files, series, mark, cc, data, label_alpha, contrast_lo, contrast_hi, color_mode, out_dir, out_file_name):
    """
    Given a list of image files, the correponding labels files, the desired series, tabular data of label features, alpha value for opacity, low contrast threshold, high contrast threshold, a color mode, the output directory as well as the output file name,
    creates an masked overlay with hover information provided for each label.
    """
    f=series-1
    image = io.imread(image_files[f])
    labels = io.imread(labels_files[f])
    s_processed_data = data.loc[data["series"] == series]
    s_processed_data = s_processed_data.reset_index(drop=True)
    fig = px.imshow(image, binary_string=True, width=800, height=800, zmin=contrast_lo, zmax=contrast_hi)
    fig.update_traces(hoverinfo='skip', hovertemplate=None)
    fig.update_layout(showlegend=False, template="plotly_dark")
    for i in tqdm(range(labels.max())):
        label_i = i+1
        contour = measure.find_contours(labels == label_i)[0] # find contour xy values for each label
        y, x = contour.T # transpose contour array to assign each x value to 'x' and each y value to 'y'
        hoverinfo = ''

        if color_mode=="cell_cycle":
            properties = ['label','series', "nuc_CTCF_n", "nuc_median_intensity_bgCorr_n", "nuc_volume", f'{mark}_CTCF_n', f'{mark}_median_intensity_bgCorr_n', f'{mark}_volume', f'{cc}_CTCF_n', f'{cc}_median_intensity_bgCorr_n', "nuc_circularity", "cell", 'cell_cycle']

            if data[(data["label"]==label_i) & (data["series"]==series)]["cell"].values[0]=="no_cell":
                c="red"
            elif data[(data["label"]==label_i) & (data["series"]==series)]["cell_cycle"].values[0]=="G1":
                c="cyan"
            elif data[(data["label"]==label_i) & (data["series"]==series)]["cell_cycle"].values[0]=="G2":
                c="magenta"

            for prop_name in properties:
                hoverinfo += f'<b>{prop_name}: {s_processed_data[prop_name][i]}</b><br>'
            fig.add_trace(go.Scatter(
                x=x, y=y, name=label_i, mode = 'lines', fill='toself', line_color=c, opacity=label_alpha, showlegend=False,
                hovertemplate=hoverinfo, hoveron='points+fills'))

        elif color_mode=="ploidy":

            properties = ['label','series', "nuc_CTCF_n", "nuc_median_intensity_bgCorr_n", "nuc_volume", f'{mark}_CTCF_n', f'{mark}_median_intensity_bgCorr_n', f'{mark}_volume', f'{cc}_CTCF_n', f'{cc}_median_intensity_bgCorr_n', "nuc_circularity", "cell", 'cell_cycle', 'ploidy']

            if data[(data["label"]==label_i) & (data["series"]==series)]["cell"].values[0]=="no_cell":
                c="red"
            elif data[(data["label"]==label_i) & (data["series"]==series)]["ploidy"].values[0]=="NA":
                c="yellow"
            elif data[(data["label"]==label_i) & (data["series"]==series)]["ploidy"].values[0]=="2N":
                c="cyan"
            elif data[(data["label"]==label_i) & (data["series"]==series)]["ploidy"].values[0]=="4N":
                c="magenta"
            for prop_name in properties:
                hoverinfo += f'<b>{prop_name}: {s_processed_data[prop_name][i]}</b><br>'
            fig.add_trace(go.Scatter(
                x=x, y=y, name=label_i, mode = 'lines', fill='toself', line_color=c, opacity=label_alpha, showlegend=False,
                hovertemplate=hoverinfo, hoveron='points+fills'))

        elif color_mode=="cell_type":

            properties = ['label','series', "nuc_CTCF_n", "nuc_median_intensity_bgCorr_n", "nuc_volume", f'{mark}_CTCF_n', f'{mark}_median_intensity_bgCorr_n', f'{mark}_volume', f'{cc}_CTCF_n', f'{cc}_median_intensity_bgCorr_n', "nuc_circularity", "cell", 'cell_cycle', 'ploidy', 'cell_type']

            
            if data[(data["label"]==label_i) & (data["series"]==series)]["cell"].values[0]=="no_cell":
                c="red"
            elif data[(data["label"]==label_i) & (data["series"]==series)]["cell_type"].values[0]==f"{mark}_positive":
                c="magenta"
            elif data[(data["label"]==label_i) & (data["series"]==series)]["cell_type"].values[0]==f"{mark}_negative":
                c="cyan"
            for prop_name in properties:
                hoverinfo += f'<b>{prop_name}: {s_processed_data[prop_name][i]}</b><br>'
            fig.add_trace(go.Scatter(
                x=x, y=y, name=label_i, mode = 'lines', fill='toself', line_color=c, opacity=label_alpha, showlegend=False,
                hovertemplate=hoverinfo, hoveron='points+fills'))
        else:

            properties = ['label','series', "nuc_CTCF_n", "nuc_median_intensity_bgCorr_n", "nuc_volume", f'{mark}_CTCF_n', f'{mark}_median_intensity_bgCorr_n', f'{mark}_volume', f'{cc}_CTCF_n', f'{cc}_median_intensity_bgCorr_n', "nuc_circularity"]

            for prop_name in properties:
                hoverinfo += f'<b>{prop_name}: {s_processed_data[prop_name][i]}</b><br>'
            fig.add_trace(go.Scatter(
                x=x, y=y, name=label_i, mode = 'lines', fill='toself', opacity=label_alpha, showlegend=False,
                hovertemplate=hoverinfo, hoveron='points+fills'))
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    fig.write_html(os.path.join(out_dir, out_file_name + f"_series{series}.html"))