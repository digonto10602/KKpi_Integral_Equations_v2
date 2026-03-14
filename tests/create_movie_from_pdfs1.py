
import os
import sys
import subprocess
from pathlib import Path

def create_movie_from_pdfs(input_prefix='plot', num_files=50, output_movie='output.mp4', 
                          fps=5, duration=10, dpi=150):
    """
    Create a movie from PDF files

    Parameters:
    -----------
    input_prefix : str
        Prefix of input PDF files (default: 'plot')
    num_files : int
        Number of PDF files to process (default: 50)
    output_movie : str
        Output movie filename (default: 'output.mp4')
    fps : int
        Frames per second (auto-calculated from duration if None)
    duration : float
        Total duration of movie in seconds (default: 10)
    dpi : int
        DPI for PDF to image conversion (default: 150)
    """

    print("="*70)
    print("PDF to Movie Converter")
    print("="*70)

    # Calculate fps from duration and number of files
    if duration:
        fps = num_files / duration
        print(f"Duration: {duration} seconds")
        print(f"Number of frames: {num_files}")
        print(f"Calculated FPS: {fps:.2f}")
    else:
        print(f"FPS: {fps}")
        print(f"Number of frames: {num_files}")
        duration = num_files / fps
        print(f"Calculated duration: {duration:.2f} seconds")

    # Check if ffmpeg is available
    try:
        subprocess.run(['ffmpeg', '-version'], stdout=subprocess.PIPE, 
                      stderr=subprocess.PIPE, check=True)
        ffmpeg_available = True
    except (subprocess.CalledProcessError, FileNotFoundError):
        ffmpeg_available = False
        print("\nWarning: ffmpeg not found. Will try alternative method.")

    # Create temporary directory for PNG files
    temp_dir = Path('temp_frames')
    temp_dir.mkdir(exist_ok=True)
    print(f"\nCreating temporary directory: {temp_dir}")

    # Convert PDFs to PNG
    print("\nConverting PDFs to PNG images...")
    missing_files = []

    try:
        from pdf2image import convert_from_path
        pdf2image_available = True
    except ImportError:
        pdf2image_available = False
        print("pdf2image not installed. Trying alternative method...")

    for i in range(1, num_files + 1):
        pdf_file = f"{input_prefix}_{i}.pdf"
        png_file = temp_dir / f"frame_{i:04d}.png"

        if not os.path.exists(pdf_file):
            missing_files.append(pdf_file)
            continue

        # Try pdf2image first
        if pdf2image_available:
            try:
                images = convert_from_path(pdf_file, dpi=dpi)
                images[0].save(png_file, 'PNG')
                print(f"  Converted: {pdf_file} -> {png_file}")
            except Exception as e:
                print(f"  Error converting {pdf_file}: {e}")
        else:
            # Try using ImageMagick's convert
            try:
                subprocess.run([
                    'convert',
                    '-density', str(dpi),
                    pdf_file,
                    '-quality', '100',
                    str(png_file)
                ], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                print(f"  Converted: {pdf_file} -> {png_file}")
            except (subprocess.CalledProcessError, FileNotFoundError) as e:
                print(f"  Error: Could not convert {pdf_file}")
                print(f"  Install either 'pdf2image' or 'ImageMagick'")

    if missing_files:
        print(f"\nWarning: {len(missing_files)} files not found:")
        for f in missing_files[:5]:
            print(f"  - {f}")
        if len(missing_files) > 5:
            print(f"  ... and {len(missing_files)-5} more")

    # Count successfully converted frames
    frame_files = sorted(temp_dir.glob('frame_*.png'))
    num_frames = len(frame_files)

    if num_frames == 0:
        print("\nError: No frames were created. Cannot make movie.")
        print("\nPlease install one of the following:")
        print("  1. pdf2image: pip install pdf2image")
        print("  2. ImageMagick: apt-get install imagemagick (Linux) or brew install imagemagick (Mac)")
        return False

    print(f"\nSuccessfully converted {num_frames} frames")

    # Create movie using ffmpeg with filters to fix dimension issues
    if ffmpeg_available:
        print(f"\nCreating movie: {output_movie}")
        try:
            # Method 1: Use scale filter to ensure even dimensions
            # The scale filter with -2 ensures dimensions are divisible by 2
            subprocess.run([
                'ffmpeg',
                '-y',  # Overwrite output file
                '-framerate', str(fps),
                '-i', str(temp_dir / 'frame_%04d.png'),
                '-vf', 'scale=trunc(iw/2)*2:trunc(ih/2)*2',  # Ensure even dimensions
                '-c:v', 'libx264',
                '-pix_fmt', 'yuv420p',
                '-crf', '23',
                output_movie
            ], check=True)
            print(f"\n✓ Movie created successfully: {output_movie}")
            print(f"  Duration: {duration:.2f} seconds")
            print(f"  FPS: {fps:.2f}")
            success = True
        except subprocess.CalledProcessError as e:
            print(f"\nError creating movie: {e}")
            print("\nTrying alternative encoding method...")

            # Method 2: Use pad filter as fallback
            try:
                subprocess.run([
                    'ffmpeg',
                    '-y',
                    '-framerate', str(fps),
                    '-i', str(temp_dir / 'frame_%04d.png'),
                    '-vf', 'pad=ceil(iw/2)*2:ceil(ih/2)*2',  # Pad to even dimensions
                    '-c:v', 'libx264',
                    '-pix_fmt', 'yuv420p',
                    '-crf', '23',
                    output_movie
                ], check=True)
                print(f"\n✓ Movie created successfully: {output_movie}")
                success = True
            except subprocess.CalledProcessError as e2:
                print(f"Error with alternative method: {e2}")
                success = False
    else:
        print("\nError: ffmpeg is required to create the movie.")
        print("Install ffmpeg:")
        print("  Ubuntu/Debian: sudo apt-get install ffmpeg")
        print("  macOS: brew install ffmpeg")
        print("  Windows: Download from https://ffmpeg.org/")
        success = False

    # Cleanup temporary files
    cleanup = input("\nDelete temporary PNG files? (y/n): ").lower()
    if cleanup == 'y':
        import shutil
        shutil.rmtree(temp_dir)
        print(f"Deleted temporary directory: {temp_dir}")
    else:
        print(f"Temporary frames saved in: {temp_dir}")

    return success


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("USAGE: python create_movie_from_pdfs.py <input_prefix> [num_files] [output_movie] [duration]")
        print()
        print("Arguments:")
        print("  input_prefix : prefix of PDF files (e.g., 'plot' for plot_1.pdf)")
        print("  num_files    : number of PDF files (default: 50)")
        print("  output_movie : output movie filename (default: 'output.mp4')")
        print("  duration     : movie duration in seconds (default: 10)")
        print()
        print("Examples:")
        print("  python create_movie_from_pdfs.py plot")
        print("  python create_movie_from_pdfs.py plot 50 my_movie.mp4")
        print("  python create_movie_from_pdfs.py plot 50 my_movie.mp4 10")
        print("  python create_movie_from_pdfs.py G_plot 30 G_animation.mp4 15")
        print()
        print("Requirements:")
        print("  - pdf2image: pip install pdf2image")
        print("    OR")
        print("  - ImageMagick: apt-get install imagemagick")
        print("  AND")
        print("  - ffmpeg: apt-get install ffmpeg")
        sys.exit(1)

    input_prefix = sys.argv[1]
    num_files = int(sys.argv[2]) if len(sys.argv) > 2 else 50
    output_movie = sys.argv[3] if len(sys.argv) > 3 else 'output.mp4'
    duration = float(sys.argv[4]) if len(sys.argv) > 4 else 10.0

    create_movie_from_pdfs(input_prefix, num_files, output_movie, duration=duration)
