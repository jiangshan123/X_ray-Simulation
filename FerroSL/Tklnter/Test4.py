import tkinter as tk

root = tk.Tk()

canvas_width = 200
canvas_height = 200
image = tk.Canvas(root, width=canvas_width, height=canvas_height)
image.pack()
image.create_rectangle(10, 150, 200, 150, fill="blue")
image.create_text(100, 175, text="substrate")
image.create_rectangle(0, 100, 200, 100, fill="red")
image.create_text(100, 125, text="electrode")
image.create_rectangle(0, 0, 200, 0, fill="green")
image.create_text(100, 50, text="superlattice")

root.mainloop()