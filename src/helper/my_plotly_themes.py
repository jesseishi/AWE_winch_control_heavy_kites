import plotly.graph_objects as go
import plotly.io as pio
from screeninfo import get_monitors

# I want figures to nicely fit on a page in my thesis without the font size getting
# distorted. So I use a fullwidth and 0.7width template that set the appropriate size.
a4_width_cm = 21
thesis_linewidth_cm = a4_width_cm * 0.75

# Plotly can only export based on pixels, so we need to convert from cm to pixels, which
# depends on the screen you're using.
for m in get_monitors():
    if m.is_primary:
        thesis_linewidth_pix = int(thesis_linewidth_cm * m.width / (m.width_mm / 10))

print(f"On this screen {thesis_linewidth_cm} cm is {thesis_linewidth_pix} pixels wide.")


# Set default image export values.
pio.kaleido.default_width = thesis_linewidth_pix
# TODO: pio.kaleido.scope.mathjax for Latex rendering of letters.


# TODO: Make a cool thesis template.
pio.templates["thesis"] = go.layout.Template(
    layout_annotations=[
        dict(
            name="draft watermark",
            text="DRAFT",
            textangle=-30,
            opacity=0.1,
            font=dict(color="black", size=100),
            xref="paper",
            yref="paper",
            x=0.5,
            y=0.5,
            showarrow=False,
        )
    ],
)

# To set as default:
# pio.templates.default = "plotly+thesis"
