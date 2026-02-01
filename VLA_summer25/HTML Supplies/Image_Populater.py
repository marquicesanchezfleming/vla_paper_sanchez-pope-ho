import os
import pandas as pd
import logging

def generate_html():
    csv_path = 'data.csv'
    df = pd.read_csv(csv_path)

    diff_days_path = "/Users/Djslime07/VLA_Summer25/HTML Supplies/Diff_days.csv"
    diff_days_df = pd.read_csv(diff_days_path, header=None)

    print(diff_days_df.head())

    diff_days_dict = {}
    for entry in diff_days_df.iloc[:, 0]:
        parts = entry.split(' has a difference in days of ')
        obj_epoch = parts[0]
        diff_days = parts[1]
        diff_days_dict[obj_epoch] = diff_days

    epochs = ["1.1v2", "1.2v2", "2.1", "2.2", "3.1", "3.2"]

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
'''

    for index, row in df.iterrows():
        transient_name = row['name']
        ra = row['ra']
        dec = row['dec']
        IAUID = row['IAUID']
        type = row['type']

        row_html = f'''
        <div class="row">
            <div class="description">
                {transient_name} <br>
                {IAUID} <br>
                Right Ascension = <br>
                {ra} <br>
                Declination = <br>
                {dec} <br>
                Type = {type} <br>
            </div>
        '''

        for epoch in epochs:
            image_filename = f"{transient_name}_VLASS{epoch}.png"
            image_path = f"new_images/{image_filename}"
            obj_epoch = f"{transient_name}_VLASS{epoch}"
            diff_days_text = diff_days_dict.get(obj_epoch, "N/A")
            if os.path.exists(image_path):
                print(f"Image found: {image_path}")
                row_html += f'''
                <div class="box">
                    <div class="box-text">VLASS Epoch {epoch}</div>
                    <img src="{image_path}" alt="{image_filename}">
                    <div class="box-text">Difference in Days: {diff_days_text}</div>
                </div>
                '''
            else:
                print(f"Image not found: {image_path}")
                

        leg_survey_image_filename = f"{transient_name}_LegSurvey.png"
        leg_survey_image_path = f"new_images/{leg_survey_image_filename}"
        if os.path.exists(leg_survey_image_path):
            print(f"LegSurvey image found: {leg_survey_image_path}")
            row_html += f'''
            <div class="box">
                <img src="{leg_survey_image_path}" alt="{leg_survey_image_filename}">
            </div>
            '''
        else:
            print(f"LegSurvey image not found: {leg_survey_image_path}")

        ps1_image_filename = f"{transient_name}_ps1.png"
        ps1_image_path = f"new_images/{ps1_image_filename}"
        if os.path.exists(ps1_image_path):
            print(f"Extra image found: {ps1_image_path}")
            row_html += f'''
            <div class="box">
                <img src="{ps1_image_path}" alt="{ps1_image_filename}">
            </div>
            '''
        else:
            print(f"PanStars image not found: {ps1_image_path}")

        row_html += '</div>'  
        html_content += row_html

    html_content += '''
    </div>
</body>
</html>
'''

    output_path = "unfortunate.html"
    with open(output_path, 'w') as file:
        file.write(html_content)

    print(f"HTML file '{output_path}' was successfully created.")

generate_html()