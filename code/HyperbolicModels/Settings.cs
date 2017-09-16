namespace HyperbolicModels
{
	using System.Collections.Generic;
	using System.Linq;
	using System.Runtime.Serialization;
	using System.Text;

	[DataContract( Namespace = "" )]
	public class Settings
	{
		public string HoneycombString
		{
			get
			{
				if( Angles.Length == 6 )
					return "Goursat domain with dihedrals " + 
						string.Join( ",", Angles );

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
				return SaveToString( Angles );
			}
			set
			{
				Angles = LoadFromString( value );
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

		internal static string SaveToString( int[] vals )
		{
			return string.Join( ",", vals );
		}

		internal static int[] LoadFromString( string str )
		{
			return str.Split( new char[] { ',' } ).Select( s => int.Parse( s ) ).ToArray();
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
		public string DisplayString
		{
			get
			{
				return string.Format(
					"Number of Edges: {0}\nEdge Width: {1}\nActiveMirrors: {2}",
					NumEdges, EdgeWidth, ActiveMirrors );
			}
		}

		/// <summary>
		/// The active mirrors, 0-indexed.  Array only contains active mirrors, so e.g.
		/// Wikipedia's 1010 would end up as { 0, 2 }
		/// </summary>
		public int[] Active { get; set; }

		/// <summary>
		/// This is only here to pretty up the xml persistence.
		/// </summary>
		[DataMember]
		private string ActiveMirrors
		{
			get
			{
				StringBuilder sb = new StringBuilder( "0000" );
				foreach( int a in Active )
					sb[a] = '1';
				return sb.ToString();
			}
			set
			{
				List<int> active = new List<int>();
				for( int i = 0; i < 4; i++ )
					if( value[i] == '1' )
						active.Add( i );
				Active = active.ToArray();
			}
		}

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

		// Other config possibilities...
		// Resolution.
		// Colors
	}
}
