namespace R3.Drawing
{
	using R3.Geometry;
	using System.Drawing;

	public class DrawUtils
	{
		static public void DrawCircle( Circle c, Graphics g, ImageSpace i )
		{
			using( Pen pen = new Pen( Color.Black, 1.0f ) )
				DrawCircle( c, g, i, pen );
		}

		static private Rectangle? Rect( Circle c, ImageSpace i )
		{
			if( double.IsInfinity( c.Radius ) )
				return null;

			Vector3D upperLeft = i.Pixel( new Vector3D( c.Center.X - c.Radius, c.Center.Y + c.Radius, 0 ) );
			double width = i.Width( c.Radius * 2 );
			double height = i.Height( c.Radius * 2 );
			Rectangle rect = new Rectangle( (int)upperLeft.X, (int)upperLeft.Y, (int)width, (int)height );
			return rect;
		}

		static public void DrawCircle( Circle c, Graphics g, ImageSpace i, Pen p )
		{
			Rectangle? rect = Rect( c, i );
			if( rect == null )
				return;
			g.DrawEllipse( p, rect.Value );
		}

		static public void DrawFilledCircle( Circle c, Graphics g, ImageSpace i, Brush b )
		{
			Rectangle? rect = Rect( c, i );
			if( rect == null )
				return;
			g.FillEllipse( b, rect.Value );
		}

		static public void DrawLine( Vector3D p1, Vector3D p2, Graphics g, ImageSpace i, Pen p )
		{
			g.DrawLine( p, VecToPoint( p1, i ), VecToPoint( p2, i ) );
		}

		static public void DrawTriangle( Mesh.Triangle triangle, Graphics g, ImageSpace i )
		{
			using( Pen pen = new Pen( Color.Black, 1.0f ) )
			{
				g.DrawLine( pen, VecToPoint( triangle.a, i ), VecToPoint( triangle.b, i ) );
				g.DrawLine( pen, VecToPoint( triangle.b, i ), VecToPoint( triangle.c, i ) );
				g.DrawLine( pen, VecToPoint( triangle.c, i ), VecToPoint( triangle.a, i ) );
			}
		}

		static private Point VecToPoint( Vector3D vec, ImageSpace i )
		{
			Vector3D temp = i.Pixel( vec );
			return new Point( (int)temp.X, (int)temp.Y );
		}
	}
}
