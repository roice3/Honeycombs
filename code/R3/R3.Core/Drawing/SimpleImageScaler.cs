using System.Drawing;
using System.Drawing.Drawing2D;

namespace R3.Core.Drawing
{
	public class SimpleImageScaler
	{
		/// <summary>
		/// Scale up an image in the simplest way (multiplying size by an integer factor).
		/// </summary>
		public static Bitmap ScaleUpBitmap( string filename, int scale )
		{
			// Load the original image
			using( Bitmap original = new Bitmap( filename ) )
			{
				int newWidth = (int)(original.Width * scale);
				int newHeight = (int)(original.Height * scale);

				// Create a new bitmap with the new dimensions
				Bitmap scaled = new Bitmap( newWidth, newHeight );

				// Use high-quality settings for resizing
				using( Graphics g = Graphics.FromImage( scaled ) )
				{
					g.InterpolationMode = InterpolationMode.NearestNeighbor;
					g.SmoothingMode = SmoothingMode.None;
					g.PixelOffsetMode = PixelOffsetMode.None;
					g.CompositingQuality = CompositingQuality.HighQuality;

					// Draw the original image, scaled to the new size
					g.DrawImage( original, 0, 0, newWidth, newHeight );
				}

				return scaled; // caller owns this Bitmap and should dispose it when done
			}
		}

		public static void ScaleUpAndSaveBitmap( string filename, int scale, string outputName )
		{
			using( Bitmap scaled = ScaleUpBitmap( filename, scale ) )
			{
				scaled.Save( outputName );
			}
		}
	}
}
