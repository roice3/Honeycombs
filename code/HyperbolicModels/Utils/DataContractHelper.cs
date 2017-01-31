namespace HyperbolicModels
{
	using System.IO;
	using System.Runtime.Serialization;
	using System.Text;
	using System.Xml;

	/// <summary>
	/// Class with useful methods for saving/loading objects.
	/// </summary>
	public class DataContractHelper
	{
		public static void SaveToXml( object obj, string filename )
		{
			using( var writer = XmlWriter.Create( filename, WriterSettings ) )
			{
				DataContractSerializer dcs = new DataContractSerializer( obj.GetType() );
				dcs.WriteObject( writer, obj );
			}
		}

		public static object LoadFromXml( System.Type objectType, string filename )
		{
			using( var reader = XmlReader.Create( filename, ReaderSettings ) )
			{
				DataContractSerializer dcs = new DataContractSerializer( objectType );
				return dcs.ReadObject( reader, verifyObjectName: false );
			}
		}

		public static string SaveToString( object obj )
		{
			StringBuilder sb = new StringBuilder();
			using( XmlWriter writer = XmlWriter.Create( sb, WriterSettings ) )
			{
				DataContractSerializer dcs = new DataContractSerializer( obj.GetType() );
				dcs.WriteObject( writer, obj );
			}
			return sb.ToString();
		}

		public static object LoadFromString( System.Type objectType, string saved )
		{
			using( StringReader sr = new StringReader( saved ) )
			using( XmlReader reader = XmlReader.Create( sr, ReaderSettings ) )
			{
				DataContractSerializer dcs = new DataContractSerializer( objectType );
				return dcs.ReadObject( reader, verifyObjectName: false );
			}
		}

		public static XmlWriterSettings WriterSettings
		{
			get
			{
				XmlWriterSettings settings = new XmlWriterSettings();
				settings.OmitXmlDeclaration = true;
				settings.Indent = true;
				return settings;
			}
		}

		public static XmlReaderSettings ReaderSettings
		{
			get
			{
				XmlReaderSettings settings = new XmlReaderSettings();
				settings.IgnoreWhitespace = true;
				return settings;
			}
		}
	}
}
