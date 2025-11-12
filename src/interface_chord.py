import random

from collections import defaultdict

from highcharts_core.chart import Chart
from highcharts_core.options import HighchartsOptions
from highcharts_core.options.series.sankey import SankeySeries

from area_interactions import load_and_preprocess_data


def simulate_residue_contacts(pivot_df, num_pairs=20, seed=42):
    totals = pivot_df.sum(axis=0).sort_values(ascending=False)
    top_interaction = totals.idxmax()

    random.seed(seed)
    contacts = defaultdict(int)

    while len(contacts) < num_pairs:
        residue_a = random.randint(1, 50)
        residue_b = random.randint(51, 100)
        pair_key = (f'Residue_{residue_a}', f'Residue_{residue_b}')
        contacts[pair_key] += random.randint(10, 100)

    contact_points = [
        {'from': pair[0], 'to': pair[1], 'weight': freq}
        for pair, freq in contacts.items()
    ]

    return top_interaction, contact_points


def generate_chord_diagram(top_interaction, contact_points, output='interface_chord_diagram.html'):
    series = [
        SankeySeries(
            name=f'{top_interaction} Residue Connectivity',
            data=contact_points,
            data_labels={'enabled': True},
            tooltip={
                'pointFormat': '<span style="color:{point.color}">{point.from}</span>'
                               ' → <span style="color:{point.color}">{point.to}</span>: '
                               '<b>{point.weight}</b>'
            },
            show_in_legend=False
        )
    ]

    options = HighchartsOptions(
        chart={
            'type': 'sankey',
            'height': 700
        },
        title={
            'text': 'Simulated Residue Interface Chord Diagram'
        },
        subtitle={
            'text': f'{top_interaction} contact frequencies between residues (Fragment_A: 1-50, Fragment_B: 51-100)'
        },
        tooltip={'useHTML': True},
        no_data={'text': 'No contacts simulated'},
        series=series
    )

    chart = Chart(options=options)
    js_literal = chart.to_js_literal()
    js_literal = js_literal.replace("Highcharts.chart(null,", "Highcharts.chart('container',")

    html_content = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="UTF-8">
    <title>Interface Chord Diagram</title>
    <script src="https://code.highcharts.com/highcharts.js"></script>
    <script src="https://code.highcharts.com/modules/sankey.js"></script>
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
        top_interaction, contact_points = simulate_residue_contacts(pivot_df)
        generate_chord_diagram(top_interaction, contact_points)
        print(f"\nTop interaction: {top_interaction}")
        print(f"Simulated residue connections: {len(contact_points)}")
    except ValueError as ve:
        print(f"\nData Loading Error: {ve}")
    except Exception as e:
        print(f"\nAn error occurred: {e}")

