from matplotlib.colors import to_rgb, to_hex

def darken_color(hex_color, factor=0.6):
    rgb = to_rgb(hex_color)
    dark_rgb = tuple(factor * c for c in rgb)
    return to_hex(dark_rgb)
