import pandas as pd
import os

def generate_html():
    csv_path = 'data.csv'
    df = pd.read_csv(csv_path)

    diff_days_path = "/Users/Djslime07/VLA_Summer25/HTML Supplies/Diff_days.csv"
    diff_days_df = pd.read_csv(diff_days_path, header=None)
    diff_days_dict = {}
    for entry in diff_days_df.iloc[:, 0]:
        parts = entry.split(' has a difference in days of ')
        obj_epoch = parts[0]
        diff_days = parts[1]
        diff_days_dict[obj_epoch] = diff_days

    band_diff_path = "/Users/Djslime07/VLA_Summer25/HTML Supplies/cutout_outputs.csv"  
    band_diff_df = pd.read_csv(band_diff_path)
    band_diff_dict = dict(zip(band_diff_df['tsn_name'], band_diff_df[' difference']))

    epochs = ["1.1v2", "1.2v2", "2.1", "2.2", "3.1", "3.2"]
    bands = ["S", "C", "X"]

    html_content = '''<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>ZTF VLASS Cross-Matching</title>
    <link rel="stylesheet" href="styles.css">
</head>
<body>
    <div class="container">
    <h1>VLASS Epochs + 25A-386 Observations and Cutouts</h1>
        <p>
            This page displays the three epochs of the VLASS plus the observations from project 25A-386. The cutouts have been normalized 
            to a common scale for visual comparison. That is, S band cutouts are 24 x 24 arcseconds, C band cutouts are 12 x 12 arcseconds, 
            and X band cutouts are 6 x 6 arcseconds. That said, it all boils down to visualization, so how it "looks" is mostly arbitrary 
            and irrelevant. The difference in days between the optical explosion date and the 25A-386 observations is also provided for each 
            transient.
        </p>
'''

    for index, row in df.iterrows():
        transient_name = row['name']
        tsn_name = row['IAUID']
        row_html = f'''
        <div class="row">
        '''
        image_count = 0  

        for epoch in epochs:
            image_filename = f"{transient_name}_VLASS{epoch}.png"
            image_path = f"new_images/{image_filename}"
            obj_epoch = f"{transient_name}_VLASS{epoch}"
            diff_days_text = diff_days_dict.get(obj_epoch, "N/A")
            if os.path.exists(image_path):
                image_count += 1
                print(f"Image found: {image_path}")
                row_html += f'''
                <div class="box">
                    <div class="box-text">{transient_name} </br> {tsn_name} </br> VLASS Epoch {epoch}</div>
                    <img src="{image_path}" alt="{image_filename}">
                    <div class="box-text"> Difference in Days: {diff_days_text}</div>
                </div>
                '''
            else:
                print(f"Image not found: {image_path}")
        
        for band in bands:
            band_filename = f"{transient_name}_{band.lower()}band.png"
            band_path = f"Images/{band_filename}"
            band_key = f"{tsn_name}_{band} Band"
            band_days = band_diff_dict.get(band_key, "N/A")

            if os.path.exists(band_path):
                image_count += 1
                row_html += f'''
                <div class="box">
                    <div class="box-text">{transient_name} </br> {tsn_name} </br> VLA {band} Band</div>
                    <img src="{band_path}" alt="{band_filename}">
                    <div class="box-text"> Difference in Days: {float(band_days):.3f}</div>
                </div>
                '''

        while image_count < 6:
            row_html += '''
            <div class="box">
                <div class="box-text">No 25A-386 Observations yet</div>
                <img src="new_images/sad_face.png" alt="placeholder image">
            </div>
            '''
            image_count += 1
        row_html += '''
        </div>
        '''
        html_content += row_html

    html_content += '''
    </div>
</body>
</html>
'''
    output_path = "success.html"
    with open(output_path, 'w') as file:
        file.write(html_content)
        print(f"HTML file '{output_path}' was successfully created.")

generate_html()
