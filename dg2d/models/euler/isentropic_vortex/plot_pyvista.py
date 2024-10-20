import pyvista as pv
pv.global_theme.font.size = 30

reader = pv.get_reader("solution.xdmf")
nt = reader.number_time_points
print("Number of solutions = ", nt)

# Last solution
reader.set_active_time_point(nt-1)
time = reader.active_time_value
print("Reading solution at time = ", time)
mesh = reader.read()

scalar = "Pressure"
contours = mesh.contour(isosurfaces=20, scalars=scalar)

p = pv.Plotter(window_size=(1200,1200))
p.add_mesh(mesh.outline(), color="k")
p.add_mesh(mesh, cmap="viridis", scalars=scalar, show_scalar_bar=False)
p.add_scalar_bar(position_x=0.2, position_y=0.0, title_font_size=35,
                 label_font_size=30)
p.add_mesh(contours, color="white", line_width=1.5)
p.add_title(scalar+", Time="+str(time))
p.view_xy()
p.camera.zoom(1.2)
p.show()
