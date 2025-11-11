import pandas as pd
from pathlib import Path
from highcharts_core.chart import Chart
from highcharts_core.options import HighchartsOptions
from highcharts_core.options.series.area import AreaSeries

BASE_DATA_PATH = Path(__file__).parent.parent / 'data'

def load_and_preprocess_data():
    all_frames_data = []
    
    for frame_id in range(1, 8):
        frame_folder = BASE_DATA_PATH / f'frame_{frame_id}'
        file_name = f'frame_{frame_id}.pd_h.pdb_A_B_summary_table.csv'
        file_path = frame_folder / file_name

        try:
            df = pd.read_csv(file_path, usecols=['Property', 'Value'])
            df['Frame'] = frame_id
            all_frames_data.append(df)
            print(f"Successfully loaded Frame {frame_id}.")
        except FileNotFoundError:
            print(f"Error: Required file not found for Frame {frame_id} at {file_path}. Skipping this frame.")
            continue
        except Exception as e:
            print(f"Error reading Frame {frame_id}: {e}")
            continue

    if not all_frames_data:
        raise ValueError("Critical Error: No data files were successfully loaded. Check file paths and existence.")

    master_df = pd.concat(all_frames_data, ignore_index=True)
    pivot_df = master_df.pivot(index='Frame', columns='Property', values='Value')
    pivot_df = pivot_df.sort_index().fillna(0)
    
    non_zero_columns = pivot_df.columns[(pivot_df != 0).any(axis=0)]
    pivot_df = pivot_df[non_zero_columns]
    
    print(f"\nFiltered to {len(non_zero_columns)} interactions with non-zero values (out of {len(master_df['Property'].unique())} total)")
    print(f"\nConsolidated Data (First 5 Rows):\n{pivot_df.head()}")
    
    return pivot_df

def create_highcharts_series(pivot_df):
    series_list = []
    
    for idx, interaction_name in enumerate(pivot_df.columns):
        data_points = pivot_df[interaction_name].reset_index().values.tolist()
        
        series = AreaSeries(
            name=interaction_name,
            data=data_points,
            xAxis=idx,
            tooltip={
                'pointFormat': (
                    '<b>{series.name}</b><br>'
                    'Frame: {point.x}<br>'
                    'Count: {point.y}'
                )
            }
        )
        series_list.append(series)
    
    return series_list

def generate_3d_area_chart(series_data, output_filename='interaction_3d_area_chart.html'):
    options = HighchartsOptions(
        chart={
            'type': 'area',
            'options3d': {
                'enabled': True,
                'alpha': 15, 
                'beta': 20, 
                'depth': 180
            },
            'height': 1000,
            'marginRight': 80 
        },
        title={
            'text': '3D Molecular Interaction Dynamics (Frames 1-7)',
            'style': {'fontSize': '22px', 'fontWeight': 'bold'}
        },
        subtitle={
            'text': f'3D Area chart showing {len(series_data)} interaction types with non-zero values over 7 simulation frames.'
        },
        xAxis=[{
            'title': {'text': 'Frame Number (X)'},
            'allowDecimals': False,
            'min': 1,
            'max': 7,
            'visible': True
        }] + [{
            'visible': False,
            'allowDecimals': False,
            'min': 1,
            'max': 7
        } for _ in range(len(series_data) - 1)],
        yAxis={
            'title': {'text': 'Interaction Count (Y)'},
            'allowDecimals': False
        },
        plotOptions={
            'area': {
                'stacking': None,
                'depth': 100,
                'marker': {'enabled': False},
                'states': {
                    'inactive': {
                        'enabled': False
                    }
                }
            }
        },
        series=series_data
    )

    chart = Chart(options=options)
    js_literal = chart.to_js_literal()
    js_literal = js_literal.replace("Highcharts.chart(null,", "Highcharts.chart('container',")
    
    html_content = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>3D Molecular Interaction Dynamics</title>
    <script src="https://code.highcharts.com/highcharts.js"></script>
    <script src="https://code.highcharts.com/highcharts-3d.js"></script>
    <script src="https://code.highcharts.com/modules/exporting.js"></script>
</head>
<body>
    <div id="container" style="width:100%; height:1000px;"></div>
    <script>
        {js_literal}
    </script>
</body>
</html>"""
    
    with open(output_filename, 'w', encoding='utf-8') as f:
        f.write(html_content)

if __name__ == '__main__':
    try:
        interaction_pivot_df = load_and_preprocess_data()
        highcharts_series_list = create_highcharts_series(interaction_pivot_df)
        generate_3d_area_chart(highcharts_series_list)
        
    except ValueError as ve:
        print(f"\nData Loading Error: {ve}")
    except Exception as e:
        print(f"\nAn error occurred during chart generation: {e}")