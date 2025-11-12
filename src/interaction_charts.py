from highcharts_core.chart import Chart
from highcharts_core.options import HighchartsOptions

from area_interactions import load_and_preprocess_data


def generate_composition_barplot(pivot_df, output='composition_barplot.html'):
    totals = pivot_df.sum(axis=0).sort_values(ascending=False)
    categories = totals.index.tolist()
    series = [{
        'name': 'Total Count',
        'data': totals.tolist()
    }]

    options = HighchartsOptions(
        chart={
            'type': 'column',
            'height': 700
        },
        title={
            'text': 'Interaction Composition Across Frames'
        },
        subtitle={
            'text': 'Total accumulated counts for all interaction types'
        },
        xAxis={
            'categories': categories,
            'title': {'text': 'Interaction Type'},
            'labels': {'rotation': -45}
        },
        yAxis={
            'title': {'text': 'Total Count'},
            'allowDecimals': False
        },
        legend={'enabled': False},
        tooltip={
            'headerFormat': '<span style="font-size: 10px">{point.key}</span><br/>',
            'pointFormat': '<span style="color:{series.color}">{series.name}</span>: <b>{point.y}</b>'
        },
        series=series
    )

    chart = Chart(options=options)
    js_literal = chart.to_js_literal()
    js_literal = js_literal.replace("Highcharts.chart(null,", "Highcharts.chart('container',")

    html_content = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>Interaction Composition Barplot</title>
    <script src="https://code.highcharts.com/highcharts.js"></script>
    <script src="https://code.highcharts.com/modules/exporting.js"></script>
</head>
<body>
    <div id="container" style="width:100%; height:700px;"></div>
    <script>
        {js_literal}
    </script>
</body>
</html>"""

    with open(output, 'w', encoding='utf-8') as f:
        f.write(html_content)

    return totals


def generate_time_traces(pivot_df, totals, output='time_traces.html'):
    top_properties = totals.head(5).index.tolist()
    frames = pivot_df.index.tolist()
    series = [{
        'name': prop,
        'data': pivot_df[prop].tolist()
    } for prop in top_properties]

    options = HighchartsOptions(
        chart={
            'type': 'spline',
            'height': 700
        },
        title={
            'text': 'Top Interactions Over Time'
        },
        subtitle={
            'text': 'Frame-by-frame counts for the five most frequent interactions'
        },
        xAxis={
            'categories': frames,
            'title': {'text': 'Frame Number'},
            'allowDecimals': False
        },
        yAxis={
            'title': {'text': 'Interaction Count'},
            'allowDecimals': False
        },
        tooltip={
            'shared': True,
            'headerFormat': '<span style="font-size: 10px">Frame {point.key}</span><br/>',
            'pointFormat': '<span style="color:{series.color}">{series.name}</span>: <b>{point.y}</b><br/>'
        },
        series=series
    )

    chart = Chart(options=options)
    js_literal = chart.to_js_literal()
    js_literal = js_literal.replace("Highcharts.chart(null,", "Highcharts.chart('container',")

    html_content = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>Interaction Time Traces</title>
    <script src="https://code.highcharts.com/highcharts.js"></script>
    <script src="https://code.highcharts.com/modules/exporting.js"></script>
</head>
<body>
    <div id="container" style="width:100%; height:700px;"></div>
    <script>
        {js_literal}
    </script>
</body>
</html>"""

    with open(output, 'w', encoding='utf-8') as f:
        f.write(html_content)


if __name__ == '__main__':
    try:
        pivot_df = load_and_preprocess_data()
        totals = generate_composition_barplot(pivot_df)
        generate_time_traces(pivot_df, totals)
    except ValueError as ve:
        print(f"\nData Loading Error: {ve}")
    except Exception as e:
        print(f"\nAn error occurred: {e}")

