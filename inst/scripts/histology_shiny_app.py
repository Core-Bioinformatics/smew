
import matplotlib.pyplot as plt
from shiny import App, ui, render, reactive, types
import numpy as np
from skimage.io import imread
import pandas as pd
import os
import glob
import copy
from skimage.transform import AffineTransform, warp
from skimage import img_as_float



## --- Read metadata and intensity files globally ---


meta_path = os.path.join('input', 'MSI_metadata.csv')
intens_path = os.path.join('input', 'MSI_intensities.csv')
if not os.path.isfile(meta_path):
    raise FileNotFoundError(f"Required metadata file not found: {meta_path}")
if not os.path.isfile(intens_path):
    raise FileNotFoundError(f"Required intensity file not found: {intens_path}")
full_meta = pd.read_csv(meta_path, index_col=0)
full_intensities = pd.read_csv(intens_path, index_col=0)

# Define the Shiny app UI

# ---- UI Layout: Main app structure ----
app_ui = ui.page_fluid(

    # App-wide CSS for centering and styling
    ui.tags.style(
        """
        .center-box {
            max-width: 800px;  /* Set the maximum width of the box */
            margin: 20px auto;    /* Center-aligns the box */
            background-color: #D1DAE2;
            padding: 15px;
            box-shadow: 0 4px 8px rgba(0,0,0,0.1); /* Optional shadow for styling */
            border-radius: 8px;
        }
        """
    ),
    

    # ---- Title and sample selection panel ----
    ui.panel_well(
        ui.layout_columns(
            ui.panel_title("SMEW: Alignment of MSI and Histology for Manual Region Annotation"),)),

    # ---- Main controls: sample, peak selection, dim reduction ----
    ui.panel_well(
        ui.h3("Sample selection:"),  
        ui.output_text("show_samples"),
        ui.input_select("pick_sample", None, choices=[],width='400px'),
        ui.output_text("show_peaks"),
        ui.h3("Pick colouring for MSI data:"),
        ui.input_select("msi_colouring", None, choices=['PC1','First 3 PCs','Individual peak'],width='400px'),
        ui.panel_conditional("input.msi_colouring=='Individual peak'",
            ui.input_selectize("peak_choice", None, choices=[],multiple=False,width='400px')),
        ui.panel_conditional("input.msi_colouring!='Individual peak'",
            ui.input_selectize("peak_choices", None, choices=[],multiple=True,width='400px')),
        ui.input_slider("point_size","Select size of point",1,300,100,width='400px'),
        ui.input_action_button("run_dimred", "Run dimensionality reduction"),
        class_="center-box"),

    # ---- Show MSI image and histology after dim reduction ----
    ui.panel_conditional("output.middleHE=='MSI histology image detected'",
                        ui.layout_columns(
                            ui.card(ui.output_plot("show_dim_red2")),
                            ui.card(ui.output_plot("show_msi_histology")),
                            row_heights="500px"
                        )),                     
    # Save dim reduction and show status
    ui.output_text("middleHE"),

    # ---- Landmark selection and registration panel ----
    ui.panel_conditional("output.middleHE=='MSI histology image detected'",
                        ui.h3("Select landmarks between MSI data and MSI histology"),
                        ui.layout_columns(
                            ui.card(ui.output_plot("plot_MSI2HE_left", click=True)),
                            ui.card(ui.output_plot("plot_MSI2HE_right", click=True)),
                            row_heights="500px"),
                        ui.layout_columns(
                            ui.card(
                                ui.panel_conditional(
                                    "input.plot_MSI2HE_left_click", 
                                    ui.output_plot("plot_MSI2HE_withselected_left", click=False))),
                            ui.card(
                                ui.panel_conditional(
                                    "input.plot_MSI2HE_right_click", 
                                    ui.output_plot("plot_MSI2HE_withselected_right", click=False))),
                            row_heights="500px"),
                        ui.h4("Recorded Landmarks"),
                        ui.output_table("coords_MSI2HE"),
                        ui.input_action_button("download_MSI2HE", "Perform coregistration and save image"),
                        ui.input_action_button("undo_MSI2HE_left_click", "Undo last MSI Image click"),
                        ui.input_action_button("undo_MSI2HE_right_click", "Undo last MSI histology click")),

    # ---- Show aligned image after registration ----
    ui.panel_conditional(
        "output.middleHE=='MSI histology image detected'",
        ui.h3("Aligned MSI histology (transformed using landmarks)"),
        ui.card(ui.output_plot("show_aligned_he"))
    )
)

# Define the Shiny app server logic

# ---- Server: All app logic and reactivity ----
def server(input, output, session):
    # Store clicked coordinates for landmarks (left: MSI, right: histology)
    clicked_coords_MSI2HE_left = reactive.Value([])
    clicked_coords_MSI2HE_right = reactive.Value([])

    # Show available samples
    # Find all valid samples that have a histology image
    def find_samples():
        # Use global metadata
        if full_meta.empty:
            return []
        samples = full_meta['Sample'].unique()
        print(samples)
        valid_samples = []
        for sample in samples:
            hist_img = glob.glob(f'input/{sample}/histology.*')
            hist_img = [x for x in hist_img if any(ext in x for ext in ['tiff','png','jpg'])]
            if hist_img:
                valid_samples.append(sample)
        valid_samples.sort()
        print(valid_samples)
        ui.update_select("pick_sample", choices=valid_samples)
        return valid_samples

    # Display the list of found samples
    @output
    @render.text
    def show_samples():
        sample_list = find_samples()
        display_string = ", ".join(sample_list)
        display_string = "Found the following samples: " + display_string
        return(display_string)

    
    # Print whether or not there is an MSI histology image (this seems to be necessary for the reactivity but could probably try to remove)
    # Show if MSI histology image is detected for selected sample
    @output
    @render.text
    @reactive.event(input.pick_sample)
    def middleHE():
        if glob.glob('input/'+input.pick_sample()+'/histology.*') != []:
            msi_histology_img = glob.glob('input/'+input.pick_sample()+'/histology.*')
            msi_histology_img = [x for x in msi_histology_img if any(ext in x for ext in ['tiff','png','jpg'])]
            if not (msi_histology_img is []):
                return('MSI histology image detected')
            else:
                return('No MSI histology image detected')
        else:
            return('No MSI histology image detected')
    
    # Update choices for dim reduction based on whether there's an MSI histology
    # Update peak selection options for current sample
    def dimred_options():
        # Use global intensity matrix, filter for sample
        sample = input.pick_sample()
        msi_intensities = full_intensities[full_meta['Sample'] == sample]
        ui.update_selectize("peak_choice", choices=list(msi_intensities.columns))
        ui.update_selectize("peak_choices", choices=list(msi_intensities.columns))
        return list(msi_intensities.columns)
    
    # Display number of peaks in selected sample
    @output
    @render.text
    def show_peaks():
        peak_list = dimred_options()
        display_string = "Found " + str(len(peak_list)) + " peaks in selected sample"
        return(display_string)

# Perform dimensionaly reduction
    # Perform dimensionality reduction and prepare MSI coordinates for plotting
    @reactive.calc
    @reactive.event(input.run_dimred)
    def msi_dimred():

        clicked_coords_MSI2HE_left.set([])
        clicked_coords_MSI2HE_right.set([])
        
        # Use global intensity and metadata, filter for sample
        if full_intensities.empty or full_meta.empty:
            return pd.DataFrame()
        sample = input.pick_sample()
        msi_intensities = full_intensities[full_meta['Sample'] == sample]
        msi_coords = full_meta[full_meta['Sample'] == sample].copy()
        min_x = np.min(msi_coords['x'])
        min_y = np.min(msi_coords['y'])
        msi_coords['x'] = [x - min_x for x in msi_coords['x']]
        msi_coords['y'] = [y - min_y for y in msi_coords['y']]
        msi_coords['x'] = [x/2 for x in msi_coords['x']]
        msi_coords['y'] = [y/2 for y in msi_coords['y']]
        if input.msi_colouring() == 'PC1':
            from sklearn.decomposition import PCA
            if list(input.peak_choices()) != []:
                msi_intensities = msi_intensities[list(input.peak_choices())]
            pca = PCA(n_components=1)
            reduction = pca.fit_transform(msi_intensities)
            msi_coords['color']=reduction
            
        if input.msi_colouring() == 'First 3 PCs':
            from sklearn.decomposition import PCA
            from sklearn.preprocessing import MinMaxScaler
            if list(input.peak_choices()) != []:
                msi_intensities = msi_intensities[list(input.peak_choices())]
            pca = PCA(n_components=3)
            reduction = pd.DataFrame(pca.fit_transform(msi_intensities))
            scaler = MinMaxScaler()
            reduction_scaled = pd.DataFrame(scaler.fit_transform(reduction), columns=reduction.columns)
            reduction_colours = reduction_scaled.values.tolist()
            def rgb_to_hex(r, g, b):
                return '#{:02x}{:02x}{:02x}'.format(r, g, b)
            reduction_colours_hex = [rgb_to_hex(int(np.round(255*x)),int(np.round(255*y)),int(np.round(255*z))) for [x,y,z] in reduction_colours]
            msi_coords['color']=reduction_colours_hex
#            ax.scatter(x=msi_coords['x'], y=msi_coords['y'], c=reduction_colours_hex,marker='.',s=input.point_size())
        if input.msi_colouring() == 'Individual peak':
            msi_coords['color'] = list(msi_intensities[input.peak_choice()])
        return(msi_coords)
    

    
    # Create a scatter plot of the MSI data after dimensionality reduction
    @reactive.calc
    @reactive.event(input.run_dimred)
    def msi_dimred_plot():
        msi_coords = msi_dimred()
        fig, ax = plt.subplots(nrows=1, ncols=1,dpi=100)  # create figure & 1 axis
        ax.margins(x=0,y=0)
        ax.scatter(x=msi_coords['x'], y=msi_coords['y'], c=msi_coords['color'],marker='.',s=input.point_size())
        fig.gca().set_aspect('equal')
        ax.set_title('MSI Image')
        ax.set_rasterization_zorder(0)
        fig.tight_layout()
        return (fig,ax)
    

    # Output plot for MSI data (for UI)
    @output
    @render.plot(height=450)
    def show_dim_red():
        # Load the images
        fig,ax = copy.deepcopy(msi_dimred_plot())
        fig.set_dpi(100)
        return fig
    
    # Output plot for MSI data (duplicate for UI layout)
    @output
    @render.plot(height=450)
    def show_dim_red2():
        # Load the images
        fig,ax = copy.deepcopy(msi_dimred_plot())
        fig.set_dpi(100)
        return fig
    
    # Show MSI histology
    # Show the MSI histology image for the selected sample
    @output
    @render.plot(height=450)
    @reactive.event(input.pick_sample)
    def show_msi_histology():
        # Create the figure and axes
        fig, ax = plt.subplots(1, 1, figsize=(10, 5))
        # Load the images
        try:
            msi_histology = imread(glob.glob('input/'+input.pick_sample()+'/histology.*')[0])
            # Display the images in two subplots
            ax.imshow(msi_histology)
            ax.set_title('MSI Histology Image')

            # Tight layout for clean display
            fig.tight_layout()
            return fig
        except :
            pass
    
     # All elements for picking landmarks between MSI and MSI histology ---------------------------

    # Show dim red 
    @output
    @render.plot(height=450)
    def plot_MSI2HE_left():
        # Plot dim
        fig,ax = copy.deepcopy(msi_dimred_plot())
        fig.set_dpi(100)
        return fig
    
    # Show MSI histology
    @output
    @render.plot(height=450)
    @reactive.event(input.pick_sample)
    def plot_MSI2HE_right():
        # Create the figure and axes
        fig, ax = plt.subplots(1, 1, figsize=(10, 5))
        # Load the images
        try:
            msi_histology = imread(glob.glob('input/'+input.pick_sample()+'/histology.*')[0])
            # Display the images in two subplots
            ax.imshow(msi_histology)
            ax.set_title('MSI Histology Image')

            # Tight layout for clean display
            fig.tight_layout()
            return fig
        except:
            pass

    # Show dim red with clicked points
    @output
    @render.plot(height=450)
    @reactive.event(input.plot_MSI2HE_left_click,input.undo_MSI2HE_left_click)
    def plot_MSI2HE_withselected_left():
        # Create the figure and axes

        fig,ax = copy.deepcopy(msi_dimred_plot())

        # Get the list of clicked coordinates
        current_coords_left = clicked_coords_MSI2HE_left()
        if current_coords_left:
            if len(current_coords_left)>0:
                x_vals, y_vals = zip(*current_coords_left)  # Unpack the coordinates
                for i in range(len(current_coords_left)):
                    ax.plot(x_vals[i], y_vals[i], 'ro', markersize=5)  # Red dots on MSI Histology image
                    ax.text(x_vals[i], y_vals[i], str(i),   color='red',fontsize=9)

        # Tight layout for clean display
        fig.set_dpi(100)
        return fig
    
    # Show MSI histology with clicked points
    @output
    @render.plot(height=450)
    @reactive.event(input.plot_MSI2HE_right_click,input.undo_MSI2HE_right_click)
    def plot_MSI2HE_withselected_right():
        # Create the figure and axes
        fig, ax = plt.subplots(1, 1, figsize=(10, 5))
        
        # Load the images
        try:
            msi_histology = imread(glob.glob('input/'+input.pick_sample()+'/histology.*')[0])

            # Display the images in two subplots
            ax.imshow(msi_histology)
            ax.set_title('MSI histology Image')

            # Get the list of clicked coordinates
            current_coords_right = clicked_coords_MSI2HE_right()
            if current_coords_right:
                if len(current_coords_right)>0:
                    x_vals, y_vals = zip(*current_coords_right)  # Unpack the coordinates
                    for i in range(len(current_coords_right)):
                        ax.plot(x_vals[i], y_vals[i], 'ro', markersize=5)  # Red dots on MSI Histology image
                        ax.text(x_vals[i], y_vals[i], str(i),color='r',fontsize=9)
            # Tight layout for clean display
            fig.tight_layout()
            return fig
        except:
            pass

    # Update clicked points 
    @reactive.Effect
    @reactive.event(input.plot_MSI2HE_left_click)
    def update_MSI2HE_leftclick():
        # Capture the click coordinates from the input
        click_info_left = input.plot_MSI2HE_left_click()
        if click_info_left is not None:
            x = click_info_left['x']  # X-coordinate of the click
            y = click_info_left['y']  # Y-coordinate of the click

            # Update the list of clicked coordinates
            current_coords_left = clicked_coords_MSI2HE_left.get()
            current_coords_left.append((x, y))  # Append the new coordinates
            clicked_coords_MSI2HE_left.set(current_coords_left)

    # Update clicked points 
    @reactive.Effect
    @reactive.event(input.plot_MSI2HE_right_click)
    def update_MSI2HE_rightclick():
        click_info_right = input.plot_MSI2HE_right_click()
        if click_info_right is not None:
            x = click_info_right['x']  # X-coordinate of the click
            y = click_info_right['y']  # Y-coordinate of the click

            # Update the list of clicked coordinates
            current_coords_right = clicked_coords_MSI2HE_right.get()
            current_coords_right.append((x, y))  # Append the new coordinates
            clicked_coords_MSI2HE_right.set(current_coords_right)

    # Undo clicked points 
    @reactive.Effect
    @reactive.event(input.undo_MSI2HE_left_click)
    def undo_MSI2HE_leftclick():
        # Capture the click coordinates from the input
        current_coords_left = clicked_coords_MSI2HE_left.get()
        current_coords_left = current_coords_left[:(len(current_coords_left)-1)]
        clicked_coords_MSI2HE_left.set(current_coords_left)

    # Undo clicked points 
    @reactive.Effect
    @reactive.event(input.undo_MSI2HE_right_click)
    def undo_MSI2HE_rightclick():
        # Capture the click coordinates from the input
        current_coords_right = clicked_coords_MSI2HE_right.get()
        current_coords_right = current_coords_right[:(len(current_coords_right)-1)]
        clicked_coords_MSI2HE_right.set(current_coords_right)

    # Make table of landmarks
    @reactive.calc
    @reactive.event(input.plot_MSI2HE_left_click,input.plot_MSI2HE_right_click,input.undo_MSI2HE_left_click,input.undo_MSI2HE_right_click)
    def coords_calc_MSI2HE():
        # Retrieve the stored coordinates
        current_coords_left = clicked_coords_MSI2HE_left.get()
        current_coords_right = clicked_coords_MSI2HE_right.get()
        
        if not current_coords_left:
            return pd.DataFrame(columns=["X_left", "Y_left", "X_right", "Y_right"])

        # Create a DataFrame to display the coordinates
        df_coords_left = pd.DataFrame(current_coords_left, columns=["X_left", "Y_left"])
        df_coords_right = pd.DataFrame(current_coords_right, columns=["X_right", "Y_right"])
        df_coords = pd.concat([df_coords_left,df_coords_right],axis=1)
        return df_coords

    # Display the clicked coordinates as a table
    @output
    @render.table
    def coords_MSI2HE():
        # Retrieve the stored coordinates
        df_coords = coords_calc_MSI2HE()
        return df_coords
    
    # Download table
    @reactive.Effect
    @reactive.event(input.download_MSI2HE)
    def download_MSI2HE():
        sample = input.pick_sample()
        output_dir = os.path.join('output', sample)
        os.makedirs(output_dir, exist_ok=True)

        ### ADDED — Compute affine transform and warp histology image
    @reactive.calc
    @reactive.event(input.download_MSI2HE)  # Run once landmarks are finalized
    def aligned_msi_histology():
        # Get landmark coordinates
        df_coords = coords_calc_MSI2HE()
        if df_coords.shape[0] < 3:
            print("Need at least 3 landmark pairs.")
            return None

        # Columns 0,1 = histology (moving image)
        # Columns 2,3 = MSI (fixed coordinate system)
        src = df_coords.iloc[:, 2:4].values
        dst = df_coords.iloc[:, :2].values

        # Load histology image
        msi_histology_path = glob.glob(f"input/{input.pick_sample()}/histology.*")[0]
        img = imread(msi_histology_path)
        img = img_as_float(img)
        # Estimate the affine transform (from histology → MSI)
        tfm = AffineTransform()
        tfm.estimate(src, dst)

        # Determine output image size based on MSI coordinate range
        x_min, x_max = dst[:, 0].min(), dst[:, 0].max()
        y_min, y_max = dst[:, 1].min(), dst[:, 1].max()

        # Output shape in pixels (height, width)
        msi_coords = msi_dimred()
        output_shape = (int(msi_coords['y'].max()), int(msi_coords['x'].max()))

        print(f"Output shape: {output_shape}")

        # Warp the image into MSI coordinate space
        warped = warp(
            img,
            inverse_map=tfm.inverse,  # Important: inverse mapping
            output_shape=output_shape,
            preserve_range=True,
            order=1  # bilinear interpolation
        )

        # Convert to uint8 for plotting/saving
        warped_uint8 = (np.clip(warped, 0, 1) * 255).astype(np.uint8)

#        downsampled_img = rescale(warped, scale=0.001, anti_aliasing=True, channel_axis=-1)

        # Save
        print('Saved!')
        return warped



        # Estimate affine transformation
#        tform = transform.estimate_transform("affine", src=src, dst=dst)

        # Apply transformation to warp the image into MSI coordinate space
#        warped = transform.warp(img, inverse_map=tform.inverse)


        return transformed_image

    ### ADDED — Display aligned image

    @output
    @render.plot(height=450)
    @reactive.event(input.download_MSI2HE)
    def show_aligned_he():
        warped = aligned_msi_histology()
        if warped is None:
            return None
        sample = input.pick_sample()
        output_dir = os.path.join('output', sample)
        os.makedirs(output_dir, exist_ok=True)
        out_path = os.path.join(output_dir, 'histology_aligned.jpg')
        plt.imsave(fname=out_path, arr=warped)
        # Show popup notification
        ui.notification_show(
            f"Aligned image saved as {out_path}",
            type="success",
            duration=4000
        )
        # Overlay MSI coordinates
        msi_coords = msi_dimred()
        fig, ax = plt.subplots(figsize=(8, 8))
        ax.imshow(warped)
        ax.scatter(msi_coords["x"], msi_coords["y"], c=msi_coords["color"], s=input.point_size(), marker=".", alpha=0.1)
        ax.set_title("Aligned histology image warped into MSI coordinate space")
        ax.set_aspect("equal")
        fig.gca().set(ylim=(0,warped.shape[0]))
        fig.tight_layout()
        return fig





# Create the Shiny app
app = App(app_ui, server)
