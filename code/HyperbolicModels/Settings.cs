namespace HyperbolicModels
{
	using System.Linq;
	using System.Runtime.Serialization;

	[DataContract( Namespace = "" )]
	public class Settings
	{
		public string HoneycombString
		{
			get
			{
				return string.Format( "{{{0},{1},{2}}}", P, Q, R );
			}
		}

		/// <summary>
		/// This is only here to pretty up the xml persistence.
		/// </summary>
		[DataMember]
		private string Dihedrals
		{
			get
			{
				return SaveToString();
			}
			set
			{
				LoadFromString( value );
			}
		}

		public int[] Angles { get; set; }

		[DataMember]
		public UhsBoundarySettings UhsBoundary { get; set; }

		[DataMember]
		public PovRaySettings PovRay { get; set; }

		public int P { get { return Angles[0]; } }
		public int Q { get { return Angles[1]; } }
		public int R { get { return Angles[2]; } }

		private string SaveToString()
		{
			return string.Join( ",", Angles );
		}

		private void LoadFromString( string str )
		{
			Angles = str.Split( new char[] { ',' } ).Select( s => int.Parse( s ) ).ToArray();
		}
	}

	[DataContract( Namespace = "" )]
	public class UhsBoundarySettings
	{
		public string DisplayString
		{
			get
			{
				return string.Format(
					"Image Width: {0}\nImage Height: {1}\nBounds: {2}", 
					ImageWidth, ImageHeight, Bounds );
			}
		}

		/// <summary>
		/// The image width, in pixels.
		/// </summary>
		[DataMember]
		public int ImageWidth { get; set; }

		/// <summary>
		/// The image height, in pixels.
		/// </summary>
		[DataMember]
		public int ImageHeight { get; set; }

		/// <summary>
		/// The bounds of the image.
		/// Default is 1.0, and larger values will effectively zoom out.
		/// </summary>
		[DataMember]
		public double Bounds { get; set; }
	}

	[DataContract( Namespace = "" )]
	public class PovRaySettings
	{
		/// <summary>
		/// The number of edges to include in the output.
		/// </summary>
		[DataMember]
		public int NumEdges { get; set; }

		/// <summary>
		/// The width of edges, in a euclidean metric at the ball origin.
		/// </summary>
		[DataMember]
		public double EdgeWidth { get; set; }


		// Other possibilities.
		// Resolution.
		// Colors
	}
}
