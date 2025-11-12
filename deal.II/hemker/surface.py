import pyvista as pv

filename = "solution.vtu"
data = pv.read(filename)

warped = data.warp_by_scalar('solution')
#warped.plot(window_size=(1200,1200), show_edges=True,
            #cmap='viridis', show_scalar_bar=False)

p = pv.Plotter(window_size=(1200,1200))
p.add_mesh(warped, show_edges=True,
           cmap='viridis', show_scalar_bar=False)
p.show()
