import { DomainItem } from './Item'

export const DomainList = props => (
  <div>
    <div>
      <button type="button" onClick={props.handleAdd}>Add</button>
    </div>
    <ul>
      {props.domains.map((domain, index) => {
        return (
          <DomainItem
            key={domain.uuid}
            data={domain}
            handleRemove={props.handleRemove(index)}
            handleChange={props.handleChange(index)}
          />
        )
      })}
    </ul>
  </div>
)
