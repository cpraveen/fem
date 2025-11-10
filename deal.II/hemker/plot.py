import pyvista as pv

filename = "solution.vtu"
data = pv.read(filename)

#data.plot(scalars="solution", cpos="xy", window_size=(1200,800),
#          show_edges=True, cmap='viridis', show_scalar_bar=True,
#          show_axes=False, zoom=1.5)

scalars = "solution"

p = pv.Plotter(window_size=(1200,800))
p.add_mesh(data, scalars=scalars, show_edges=True, show_scalar_bar=False)
p.add_scalar_bar(title=scalars,
                 position_x=0.2, position_y=0.05,
                 title_font_size=35, label_font_size=35)
p.view_xy()
p.camera.zoom(1.5)

p.show()
