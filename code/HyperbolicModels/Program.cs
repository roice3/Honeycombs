namespace HyperbolicModels
{
	using System.Collections.Generic;
	using System.IO;
	using R3.Core;
	using R3.Geometry;

	class Program
	{
		static void Main( string[] args )
		{
			try
			{
				bool wiki = false;
				if( wiki )
				{
					//HoneycombGen.GoursatSet();
					//HoneycombGen.ParacompactSet();

					ViewPath path = new ViewPath();
					Vector3D[] locations = new Vector3D[]
					{
						new Vector3D(0,.1,-.5),
						new Vector3D(-.4,0,-.1),
						new Vector3D(0,-.2,-.5),
						new Vector3D(.1,0,-.8),
						new Vector3D(0,.1,-.5),
						new Vector3D(-.4,0,-.1),
						new Vector3D(0,-.2,-.5),
						new Vector3D(.1,0,-.8),
						new Vector3D(0,.1,-.5),
					};
					List<Vector3D> lookAts = new List<Vector3D>();
					path.Initialize( locations, lookAts );
					HoneycombGen.ViewPath = path;

					for( double t = 0; t < 1; t += .0016 )
					{
						//HoneycombGen.ParacompactAnimationFrame( t );
						path.Step++;
					}

					// Other settings to make configurable...
					// width
					// resolution
					// 

					HoneycombGen.ViewPath = null;
					HoneycombDef def = new HoneycombDef() { P = 3, Q = 4, R = 4 };
					int[] active = new int[] { 1 };
					int baseHue = 220;
					HoneycombGen.Paracompact( def, active, baseHue );
				}
				else
				{ 
					HoneycombPaper.DoStuff( args );
				}
			}
			catch( System.Exception ex )
			{
				System.Diagnostics.Trace.WriteLine( ex.Message + "\n" + ex.StackTrace );
				System.Console.WriteLine( ex.Message + "\n" + ex.StackTrace );
			}
		}
	}
}
