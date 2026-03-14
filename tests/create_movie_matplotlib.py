
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from pdf2image import convert_from_path
from PIL import Image
import sys

def create_movie_matplotlib(input_prefix='plot', num_files=50, 
                           output_movie='output.mp4', duration=10, dpi=150):
    """
    Create movie using matplotlib animation (alternative method)

    Parameters:
    -----------
    input_prefix : str
        Prefix of input PDF files
    num_files : int
        Number of PDF files
    output_movie : str
        Output movie filename
    duration : float
        Duration in seconds
    dpi : int
        DPI for PDF conversion
    """

    print("Loading PDF files and converting to images...")

    frames = []
    for i in range(1, num_files + 1):
        pdf_file = f"{input_prefix}_{i}.pdf"
        try:
            # Convert PDF to image
            images = convert_from_path(pdf_file, dpi=dpi)
            frames.append(np.array(images[0]))
            print(f"  Loaded: {pdf_file}")
        except Exception as e:
            print(f"  Error loading {pdf_file}: {e}")

    if not frames:
        print("\nError: No frames loaded. Cannot create movie.")
        return

    print(f"\nLoaded {len(frames)} frames")

    # Calculate interval (milliseconds per frame)
    interval = (duration * 1000) / len(frames)
    fps = len(frames) / duration

    print(f"Creating animation...")
    print(f"  Duration: {duration} seconds")
    print(f"  FPS: {fps:.2f}")
    print(f"  Interval: {interval:.2f} ms/frame")

    # Create figure
    fig, ax = plt.subplots(figsize=(12, 9))
    ax.axis('off')

    # Initialize with first frame
    im = ax.imshow(frames[0])

    def update(frame_num):
        im.set_array(frames[frame_num])
        return [im]

    # Create animation
    anim = animation.FuncAnimation(fig, update, frames=len(frames),
                                  interval=interval, blit=True, repeat=True)

    # Save animation
    print(f"\nSaving movie to: {output_movie}")
    try:
        Writer = animation.writers['ffmpeg']
        writer = Writer(fps=fps, bitrate=1800)
        anim.save(output_movie, writer=writer, dpi=100)
        print(f"\n✓ Movie created successfully: {output_movie}")
    except Exception as e:
        print(f"\nError saving movie: {e}")
        print("Make sure ffmpeg is installed.")

    plt.close()


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("USAGE: python create_movie_matplotlib.py <input_prefix> [num_files] [output_movie] [duration]")
        print()
        print("Example:")
        print("  python create_movie_matplotlib.py plot 50 output.mp4 10")
        sys.exit(1)

    input_prefix = sys.argv[1]
    num_files = int(sys.argv[2]) if len(sys.argv) > 2 else 50
    output_movie = sys.argv[3] if len(sys.argv) > 3 else 'output.mp4'
    duration = float(sys.argv[4]) if len(sys.argv) > 4 else 10.0

    create_movie_matplotlib(input_prefix, num_files, output_movie, duration, dpi=150)
