﻿namespace HyperbolicModels
{
	using System.Drawing;
	using System.Drawing.Imaging;
	using System.IO;

	/// <summary>
	/// This will take a list of images and resize/combine them into a single image grid.
	/// </summary>
	public class ImageGrid
	{
		public class Settings
		{
			public Settings()
			{
				Columns = Rows = 5;
				Width = Height = 3000;
				FileName = "output.png";
			}

			public int Columns { get; set; }
			public int Rows { get; set; }
			public int Width { get; set; }
			public int Height { get; set; }
			public string Directory { get; set; }

			/// <summary>
			/// The filenames of the input images (just leaf names).
			/// We will fill out result by rows.
			/// </summary>
			public string[] InputImages { get; set; }

			/// <summary>
			/// The output file.
			/// </summary>
			public string FileName { get; set; }
		}

		public void Generate( Settings s )
		{
			Bitmap image = new Bitmap( s.Width, s.Height );
			Graphics g = Graphics.FromImage( image );
			g.Clear( Color.White );

			int tileWidth = s.Width / s.Columns;
			int tileHeight = s.Height / s.Columns;
			Size tileSize = new Size( tileWidth, tileHeight );

			int currentRow = 0, currentCol = 0;
			foreach( string imageName in s.InputImages )
			{
				string fullFileName = Path.Combine( s.Directory, imageName );
				Bitmap tile = new Bitmap( fullFileName );

				// Resize
				tile = new Bitmap( tile, tileSize );

				// Copy to location.
				for( int i=0; i<tile.Width; i++ )
				for( int j=0; j<tile.Height; j++ )
				{
					Color c = tile.GetPixel( i, j );
					image.SetPixel( currentCol * tileWidth + i, currentRow * tileHeight + j, c );
				}

				currentCol++;
				if( currentCol >= s.Columns )
				{
					currentCol = 0;
					currentRow++;
				}
				if( currentRow >= s.Rows )
					break;
			}

			image.Save( s.FileName, ImageFormat.Png );
		}
	}
}
