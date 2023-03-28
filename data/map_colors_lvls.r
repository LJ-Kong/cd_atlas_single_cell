## mapping detailed anno2 to higher levels of type
cell_types_anno <- read.delim("/color_panel_cell_manu_groups.txt")


## mapping between main celltypes to detailed celltypes
cell_colors_overall <- c("T cells"="#fa117a", "B cells"="#ff8f8f", "Myeloid cells"="#51c4cf",  "Stromal cells"="#5149b2", "Epithelial cells"="#ffeca1", "Cycling cells"="gray")

## other mappings
disease_colors <- c(Heal = "#95d962", NonI = "#628cd9", Infl = "#d9628c")
compart_colors <- c("Immune cells" = "#ff7e0d", "Stromal cells" = "#8c564b", "Epithelial cells" = "#c5b1d5")

# for donors
godsnot_64 = c(
    "#A1C299",
    "#300018",
    "#0AA6D8",
    "#013349",
    "#00846F",
    "#372101",
    "#FFB500",
    "#C2FFED",
    "#A079BF",
    "#CC0744",
    "#C0B9B2",
    "#C2FF99",
    "#001E09",
    "#00489C",
    "#6F0062",
    "#0CBD66",
    "#EEC3FF",
    "#456D75",
    "#B77B68",
    "#7A87A1",
    "#788D66",
    "#885578",
    "#FAD09F",
    "#FF8A9A",
    "#D157A0",
    "#BEC459",
    "#456648",
    "#0086ED",
    "#886F4C",
    "#34362D",
    "#B4A8BD",
    "#00A6AA",
    "#452C2C",
    "#636375",
    "#A3C8C9",
    "#FF913F",
    "#938A81",
    "#575329",
    "#00FECF",
    "#B05B6F",
    "#8CD0FF",
    "#3B9700",
    "#04F757",
    "#C8A1A1",
    "#1E6E00",
    "#7900D7",
    "#A77500",
    "#6367A9",
    "#A05837",
    "#6B002C",
    "#772600",
    "#D790FF",
    "#9B9700",
    "#549E79",
    "#FFF69F",
    "#201625",
    "#72418F",
    "#BC23FF",
    "#99ADC0",
    "#3A2465",
    "#922329",
    "#5B4534",
    "#FDE8DC",
    "#404E55",
    "#0089A3",
    "#CB7E98",
    "#A4E804",
    "#324E72",
    "#6A3A4C"
)
